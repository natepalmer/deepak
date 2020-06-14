import itertools
import os

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from deepak.sfmap import sfmap_plot


def replace_wt(df, seq, values):
    start = df.first_valid_index()
    for i, aa in enumerate(seq):
        df.loc[i + start, aa] = values
    return df


def correlation_plot(df1, df2, sample_name, min_counts):
    x = df1["counts"]
    y = df2["counts"]
    rsq = np.corrcoef(x, y)[0, 1]
    fig, axis = plt.subplots(1, 2, figsize=(10, 5))
    ax = sns.regplot(np.log(x), np.log(y), fit_reg=False, scatter_kws={"alpha": 0.2}, ax=axis[0])
    ax.set_xlabel("Replicate 1")
    ax.set_ylabel("Replicate 2")
    ax.set_title("log read counts")
    ax.text(0.1, 0.9, r"$R^2$ = " + str(round(rsq ** 2, 3)), transform=ax.transAxes)
    a = (df1["edited_counts"] / x).mask(x < min_counts)
    b = (df2["edited_counts"] / y).mask(y < min_counts)
    select = b.notnull() & a.notnull()
    rsq2 = np.corrcoef(a[select], b[select])[0, 1]
    ax2 = sns.regplot(a, b, fit_reg=False)
    ax2.text(0.8, 0.9, r"$R^2$ = " + str(round(rsq2 ** 2, 3)), transform=ax2.transAxes)
    ax2.set_xlabel("Replicate 1")
    ax2.set_ylabel("Replicate 2")
    ax2.set_title("editing rates")
    fig.suptitle("{} replicate correlation".format(sample_name))
    plt.axis('square')
    plt.tight_layout()
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


def hmap_plot(file, df, style, wt_seq, **kwargs):
    pp = PdfPages(file)
    try:
        sfmap_plot(df=df, pdf=pp, style=style, wt=wt_seq, dimensions='wide', **kwargs)
    finally:
        pp.close()
    return


def make_heatmaps(sample, density, geom, log2_fold_change, z_scores, std_err, wt_aa_seq, min_counts, fig_dir):
    hmap_plot(fig_dir+"/{}_density.pdf".format(sample), density, "logcounts", wt_aa_seq,
              title="{} count density".format(sample))
    hmap_plot(fig_dir+"/{}_geom_rates.pdf".format(sample), geom.mask(density < min_counts), "scores", wt_aa_seq,
              title="{} log2 fold change in editing rate (geometric mean)".format(sample), vmin=-4)
    hmap_plot(fig_dir+"/{}_rates.pdf".format(sample), log2_fold_change.mask(density < min_counts), "scores", wt_aa_seq,
              title="{} log2 fold change in editing rate (arithmetic mean)".format(sample), vmin=-4)
    hmap_plot(fig_dir+"/{}_z-scores.pdf".format(sample), z_scores, "scores", wt_aa_seq,
              title="{} Z-scores with standard error".format(sample), vmin=-4, vmax=10, df_se=std_err)
    return


def make_fig_dir(sample, base_dir, append):
    if not os.path.isdir(base_dir+"figures"):
        os.mkdir(base_dir+"figures")
    fig_dir = base_dir+"figures/"+sample+append
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