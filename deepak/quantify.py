import os

import dill
import numpy as np
import pandas as pd
from Bio import SeqIO
import scipy.stats as st

import deepak.globals
import deepak.utilities
from deepak.plot import replace_wt, all_correlations, make_heatmaps, make_fig_dir

pad = 948
target_T3 = ":917*ag"
target_G3 = ":932*ag"
target_T5 = ":50*ag"
target_G5 = ":41*ag"


# MAYBE: calculate common mutations, place in separate data structure


class Quantification:

    def __init__(self, config_file):
        self.config_file = config_file

    def configure(self, config_file):
        with open(config_file) as config:
            for line in config:
                attribute, value = line.split()

    def generate_library(self):
        self.PafParser.generate_library(self.library_fn, self.alignment_position)
        return True


def load_pickled_data(fn):
    with open(fn, mode="rb") as infile:
        analysis = dill.load(infile)
    return analysis


def search_snp_paf(paf_record, target):
    target_fields = deepak.utilities.chunk_paf(target)
    assert len(target_fields) == 2  # Should be ":n*{ref}{var}"
    target_loc = int(target_fields[0][1:])
    location = 0
    fields = deepak.utilities.chunk_paf(paf_record)
    for i, item in enumerate(fields):
        if location == target_loc and item == target_fields[1]:
            return True
        elif item[0] == ":":
            location += int(item[1:])
        elif item[0] == "*":
            location += 1
        else:
            raise Exception("Disallowed character in CS string, could be indel")
    return False


def detect_edits(item, target):
    """
    Count reads in item which contain target mutation.
    Returns the number of reads containing the target mutation and the total number of reads

    *item* is a length 2 tuple comprising a library entry in the form (*name*, *list of PafRecords or cs strings*)
    *target* is a cs string specifying the target mutation to search for
    """
    name = item[0]
    edits = list(map(search_snp_paf, item[1], [target] * len(item[1])))
    return np.sum(edits), len(edits)


def decode_paf(paf_str):
    global reference_fn, pad
    ref = SeqIO.read(reference_fn, "fasta")
    fields = deepak.utilities.chunk_paf(paf_str)
    dna_loc = int(fields[0][1:])
    pos = (dna_loc + pad) // 3
    result_dna = ref[:dna_loc]
    for mut in fields[1:]:
        if mut.startswith("*"):
            result_dna += mut[2]
            dna_loc += 1
        else:
            n = int(mut[1:])
            result_dna += ref[dna_loc:dna_loc + n]
            dna_loc += n
    if dna_loc < len(ref):
        result_dna += ref[dna_loc:]
    aa = result_dna.translate()[pos - (pad // 3)]
    return int(pos), aa


def add_seq_info(data_frame):
    positions, amino_acids = list(zip(*map(decode_paf, data_frame["name"])))
    data_frame["position"] = positions
    data_frame["amino_acid"] = amino_acids
    return data_frame


def read_analysis(analysis_obj, target_mutation):
    data = {"name": [], "edited_counts": [], "counts": []}
    for member in analysis_obj.library.items():
        edited, counts = detect_edits(member, target_mutation)
        data["name"].append(member[0])
        data["edited_counts"].append(edited)
        data["counts"].append(counts)

    df = pd.DataFrame(data)
    wt = df.loc[df.name == "wt"]
    df = df.loc[df.name != "wt"]
    return df, wt


def z(p, n, wt_rate, wt_n, pooled=True, size=1):
    if n < size:
        return np.nan
    if pooled:
        combined_p = (wt_rate * wt_n + n * p) / (n + wt_n)
        return (p - wt_rate) / np.sqrt(combined_p * (1 - combined_p) * ((1 / n) + (1 / wt_n)))
    return (p - wt_rate) / np.sqrt((wt_rate * (1 - wt_rate) / wt_n) + (p * (1 - p) / n))


def add_stats(df, wt_rate, wt_n):
    n = 0
    while True:
        if "rep"+str(n)+"_counts" not in df.columns:
            break
        n += 1
    x_bar = 1
    for i in range(1, n):
        rep = "rep"+str(i)+"_"
        # Zero total counts results in NaN
        p = df[rep+"counts"]/df["counts"]
        # Members with zero counts in one replicate default to rate of other replicate, i.e. NaN ** 0 == 1
        r = (df[rep+"edited_counts"]/df[rep+"counts"]).fillna(0)
        x_bar *= np.power(r, p)
    df["geom_editing_rate"] = x_bar
    df["editing_rate"] = df["edited_counts"] / df["counts"]
    df["z-score"] = list(map(z, df["editing_rate"], df["counts"], [wt_rate] * len(df.index), [wt_n] * len(df.index)))
    df["p-value"] = st.norm.sf(np.abs(df["z-score"])) * 2  # two-tailed test
    combined_p = (wt_rate * wt_n + df["editing_rate"] * df["counts"]) / (df["counts"] + wt_n)
    df["std_error"] = np.sqrt(combined_p * (1 - combined_p) * ((1 / df["counts"]) + (1 / wt_n)))
    return df


def reference_aa(df, reference):
    start = df["position"].min()
    end = df["position"].max()
    ref = SeqIO.read(reference, "fasta")
    wt_aa_seq = str(ref.translate()[int(start - pad // 3):int(end - pad // 3) + 1].seq)
    return wt_aa_seq


def fill_aa_seq(df_seq, wt_aa_seq):
    x = set(df_seq["position"])
    least = min(x)
    y = set(range(least, least+len(wt_aa_seq)))
    z = x.difference(y)
    while len(z) > 0:
        item = z.pop()
        new_row = pd.DataFrame({"position": [item]*20, "amino_acid": deepak.globals.AA_LIST})
        df_seq = pd.concat([df_seq, new_row], sort=False, ignore_index=True)
    return df_seq


def get_plotting_frame(df, values):
    x = df.pivot(index="position", columns="amino_acid", values=values)
    if values == "counts":
        x.fillna(value=0)
    return x[deepak.globals.AA_LIST]


def read_data_set(obj_fn, target):
    analysis = load_pickled_data(obj_fn)
    df, wt = read_analysis(analysis, target)
    print("Loaded pickled data from {}".format(obj_fn))
    return df, wt


def aggregate_data(base, sample, n_replicates, target, append=""):
    obj_fns = [base.format(sample+str(i)+append) for i in range(1, n_replicates+1)]
    data_sets = list(map(read_data_set, obj_fns, [target]*n_replicates))
    return data_sets


def combine_replicates(data_sets):
    df = data_sets[0][0].copy()
    wt = data_sets[0][1].copy()
    for i, d in enumerate(data_sets):
        if i >= 1:
            # Keeps columns that are in df but not d[0] unlike df += d[0]
            df = df.combine(d[0], lambda x, y: x+y if np.issubdtype(x.dtype, np.number) else x, overwrite=False)
            wt += d[1]  # Does not have any disjoint columns
        for new_col in ("edited_counts", "counts"):
            name = "rep"+str(i)+"_"+new_col
            df[name] = d[0][new_col]
            wt[name] = d[1][new_col]
    return df, wt


def load_replicate_data(sample, base_dir, n_reps, reference, append):
    global reference_fn
    reference_fn = reference
    if not base_dir.endswith("/"):
        base_dir += "/"
    base = base_dir+"Pafparser-{}_aln/workspace.pyobj"

    if "T" in sample:
        if "3" in sample:
            target = target_T3
        else:  # 5
            target = target_T5
    else:  # G
        if "3" in sample:
            target = target_G3
        else:  # 5
            target = target_G5

    # Load data (2 replicates)
    data_sets = aggregate_data(base, sample, n_reps, target, append=append)

    df, wt = combine_replicates(data_sets)

    return df, wt, data_sets


def calculate(df, wt, reference):
    df = add_seq_info(df)

    wt_n = int(wt["counts"])
    wt_rate = float(wt["edited_counts"] / wt_n)
    wt_aa_seq = reference_aa(df, reference)

    df = add_stats(df, wt_rate, wt_n)

    density = replace_wt(get_plotting_frame(df, "counts"), wt_aa_seq, wt_n)

    log_pad = 0.000001
    geom = get_plotting_frame(df, "geom_editing_rate")
    geom_norm = replace_wt(geom / wt_rate, wt_aa_seq, 1)
    geom_fold_change = np.log2(geom_norm + log_pad)

    rates = get_plotting_frame(df, "editing_rate")
    normalized_rates = replace_wt(rates / wt_rate, wt_aa_seq, 1)
    log2_fold_change = np.log2(normalized_rates + log_pad)

    z_scores = replace_wt(get_plotting_frame(df, "z-score"), wt_aa_seq, 0)
    std_err = replace_wt(get_plotting_frame(df, "std_error"), wt_aa_seq, np.nan)

    return wt_aa_seq, density, geom_fold_change, log2_fold_change, z_scores, std_err


def pafparser_to_csv(sample, base_dir, n_reps, reference, append):
    df, wt, data_sets = load_replicate_data(sample, base_dir, n_reps, reference, append)
    for i in range(2):
        for item in ("counts", "edited_counts"):
            wt["rep{}_{}".format(i, item)] = data_sets[i][1][item]
    df_seq = add_seq_info(df)
    wt_aa_seq = reference_aa(df, reference)
    df_seq = fill_aa_seq(df_seq, wt_aa_seq)
    full = pd.concat([df_seq, wt], sort=False, ignore_index=True)
    full.to_csv(base_dir+sample+".csv")
    return df_seq, wt


def csv_to_df_wt(fn):
    df = pd.read_csv(fn, header=0, index_col=0)
    wt = df.loc[df["name"] == "wt"]
    df = df.loc[df["name"] != "wt"]
    return df, wt


def run_from_pickle(sample, base_dir, n_reps, reference, append, min_counts=1):
    if os.path.isfile(base_dir+sample+".csv"):
        df, wt = csv_to_df_wt(base_dir+sample+".csv")
    else:
        df, wt = pafparser_to_csv(sample, base_dir, n_reps, reference, append)
    wt_aa_seq, density, geom_fold_change, log2_fc, z_scores, std_err = calculate(df, wt, reference)
    fig_dir = make_fig_dir(sample, base_dir, append)
    all_correlations(df, sample, fig_dir, min_counts)
    make_heatmaps(sample, density, geom_fold_change, log2_fc, z_scores, std_err, wt_aa_seq, min_counts, fig_dir)
    return df, wt


def run_from_csv(filenames, lib_file, target_file):
    #df_0 = read_csv_to_lib_df(filenames[0], 0, library, target)
    #for i, fn in enumerate(filenames):
    pass


def read_csv_to_lib_df(fn, rep_number, library, target):
    read_df = pd.read_csv(fn, index_col=0)
    data = list()
    for member in library:
        item = (member, read_df.loc[read_df.lib_identity == member]["cs_tag"])
        edits, counts = detect_edits(item, target)
        row = {"name": member, "rep{}_edited_counts".format(rep_number): edits,
               "rep{}_counts".format(rep_number): counts}
        data.append(row)
    return pd.DataFrame(data)
