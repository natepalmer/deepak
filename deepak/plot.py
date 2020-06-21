import itertools
import os

import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from deepak.sfmap import sfmap_plot


def replace_wt(df, seq, values):
    start = df.first_valid_index()
    for i, aa in enumerate(seq):
        df.loc[i + start, aa] = values
    return df


def flatten_to_series(df):
    return pd.Series(df.to_numpy().flatten())


def correlation_plot(quant_1, quant_2, min_counts):
    x = flatten_to_series(quant_1.counts)
    y = flatten_to_series(quant_2.counts)
    rsq = np.corrcoef(x, y)[0, 1]
    fig, axis = plt.subplots(1, 2, figsize=(10, 5))
    ax = sns.regplot(np.log(x), np.log(y), fit_reg=False, scatter_kws={"alpha": 0.2}, ax=axis[0])
    ax.set_xlabel(quant_1.name)
    ax.set_ylabel(quant_2.name)
    ax.set_title("log read counts")
    ax.text(0.1, 0.9, r"$R^2$ = " + str(round(rsq ** 2, 3)), transform=ax.transAxes)
    a = (flatten_to_series(quant_1.edits) / x).mask(x < min_counts)
    b = (flatten_to_series(quant_2.edits) / y).mask(y < min_counts)
    select = b.notnull() & a.notnull()
    rsq2 = np.corrcoef(a[select], b[select])[0, 1]
    ax2 = sns.regplot(a, b, fit_reg=False)
    ax2.text(0.1, 0.9, r"$R^2$ = " + str(round(rsq2 ** 2, 3)), transform=ax2.transAxes)
    ax2.set_xlabel(quant_1.name)
    ax2.set_ylabel(quant_2.name)
    ax2.set_title("editing rates")
    fig.suptitle(f"{quant_1.name} vs {quant_2.name} replicate correlation")
    plt.axis('square')
    ax2.set_xlim(0, 0.5)
    ax2.set_ylim(0, 0.5)
    #plt.tight_layout()
    return fig


def get_sub_replicate(df, i):
    d1 = df[["rep{}_counts".format(i), "rep{}_edited_counts".format(i)]]
    d1.columns = ["counts", "edited_counts"]
    return d1


def all_correlations(df, sample_name, fig_dir, min_counts):
    counter = 0
    while "rep{}_counts".format(counter) in df:
        counter += 1
    for i, j in itertools.combinations(range(counter), 2):
        frames = [get_sub_replicate(df, x) for x in (i, j)]
        corr = correlation_plot(frames[0], frames[1], sample_name, min_counts)
        corr.savefig(fig_dir+"/{}-{}_correlation.pdf".format(sample_name, str(counter)))
        counter += 1
    return


def make_correlation_plots(quant_list, fig_dir, min_counts=1):
    for x, y in itertools.combinations(quant_list, 2):
        figure = correlation_plot(x, y, min_counts=min_counts)
        figure.savefig(fig_dir+f'/{x.name}_{y.name}_correlation.pdf')
    return


def hmap_plot(file, df, style, wt_seq, **kwargs):
    pp = PdfPages(file)
    try:
        sfmap_plot(df=df, pdf=pp, style=style, wt=wt_seq, dimensions='wide', **kwargs)
    finally:
        pp.close()
    return


def make_heatmaps(sample, stat_dict, wt_aa_seq, fig_dir, min_counts=1, lfc=True):
    start_idx = stat_dict["mean"].index[0]
    wt = wt_aa_seq[start_idx:]
    hmap_plot(fig_dir+"/{}_density.pdf".format(sample), stat_dict["counts"], "logcounts", wt,
              title="{} count density".format(sample))
    if lfc:
        pad = 0.00001
        log2_fold_change = np.log2((stat_dict["mean"]+pad) / stat_dict["mean"].loc[start_idx, wt_aa_seq[start_idx]])
        hmap_plot(fig_dir+"/{}_lfc.pdf".format(sample), log2_fold_change.mask(stat_dict["counts"] < min_counts),
                  "scores", wt, title="{} mean log2 fold change in editing rate".format(sample), vmin=-4)

        geom_lfc = np.log2((stat_dict["geom"]+pad) / stat_dict["geom"].loc[start_idx, wt_aa_seq[start_idx]])
        hmap_plot(fig_dir+"/{}_geom_lfc.pdf".format(sample), geom_lfc.mask(stat_dict["counts"] < min_counts), "scores",
                  wt, title="{} geometric mean log2 fold change in editing rate".format(sample), vmin=-4)

    else:
        # Need to create custom color map for plotting pure rates
        raise Exception("Pure rate plotting not implemented")

    hmap_plot(fig_dir+"/{}_z-scores.pdf".format(sample), stat_dict["z-scores"], "scores", wt,
              title="{} Z-scores with standard error".format(sample), vmin=-4, vmax=10, df_se=stat_dict["std_error"])
    return


def make_fig_dir(sample, base_dir, append=""):
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)
    fig_dir = os.path.join(base_dir, "figures")
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
    else:
        print("Overwriting files in {}".format(fig_dir))
    return fig_dir


def run(sample, base_dir, n_reps, reference, append, min_counts=1, target=None, basename="Pafparser-"):
    #if os.path.isfile(base_dir+sample+".csv"):
    #    df, wt = csv_to_df_wt(base_dir+sample+".csv")
    #else:
    df, wt = pafparser_to_csv(sample, base_dir, n_reps, reference, append, target=target, basename=basename)
    wt_aa_seq, density, geom_fold_change, log2_fc, z_scores, std_err = calculate(df, wt, reference)
    fig_dir = make_fig_dir(sample, base_dir, append)
    all_correlations(df, sample, fig_dir, min_counts)
    make_heatmaps(sample, density, geom_fold_change, log2_fc, z_scores, std_err, wt_aa_seq, min_counts, fig_dir)
    return df, wt