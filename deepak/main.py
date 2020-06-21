import os
import argparse

from deepak.experiment import Experiment
from deepak.utilities import make_output_directory
from deepak.quantify import run_files


def analyze(**kwargs):
    paf_fn = kwargs["aligned_reads"]
    library_fn = kwargs["library"]
    reference_fn = kwargs["reference"]
    alignment_pos = kwargs["position"]
    if kwargs["output"]:
        out_dir = make_output_directory(kwargs["output"])
    else:
        out_dir = make_output_directory(paf_fn)

    analysis = Experiment(paf_fn, reference_fn)
    analysis.library.construct(library_fn, alignment_pos)
    if kwargs["excluded"] is not None:
        with open(kwargs["excluded"]) as mut_file:
            excluded = [x.strip() for x in mut_file.readlines()]
        analysis.set_excluded_mutations(excluded)
    if kwargs["filter"] is not None:
        analysis.set_filters({"substitution": int(kwargs["filter"])})
    analysis.run(out_dir)
    #analysis.export(out_dir)
    return analysis


def build_classify_parser():
    parser = argparse.ArgumentParser(description="Classify amplicon library data for deep mutational scan")
    parser.add_argument("reads", help="Reads aligned in .paf format generated by minimap2 with the --cs flag")
    parser.add_argument("-l", "--library", required=True,
                        help="File with each line corresponding to the sequence of a library element")
    parser.add_argument("-r", "--reference", required=True, help="Reference sequence used for alignment")
    parser.add_argument("-p", "--position", default=0, type=int,
                        help="Location of start of library sequences within reference")
    parser.add_argument("-o", "--output", default=False, help="Directory to write output files")
    parser.add_argument("-e", "--excluded", default=None,
                        help="File with each line containing the cs string of a mutation to be excluded from filtering")
    parser.add_argument("-f", "--filter", default=0, help="Max number of non-library mutations allowed")
    return parser


def build_quantify_parser():
    parser = argparse.ArgumentParser(description="Analyze amplicon library data for deep mutational scan")
    parser.add_argument("files", nargs="+",
                        help="csv files containing the valid reads generated by deepak classify, one file per replicate")
    parser.add_argument("-n", "--name", required=True, help="Name of the sample to be used in creating output files")
    parser.add_argument("-l", "--library", required=True,
                        help="File with each line corresponding to the sequence of a library element")
    parser.add_argument("-r", "--reference", required=True, help="Reference sequence used for alignment")
    parser.add_argument("-p", "--position", default=0, type=int,
                        help="Location of start of library sequences within reference")
    parser.add_argument("-o", "--output", required=True, help="Directory to write output files")
    parser.add_argument("-t", "--target", required=True, help="Target mutation in paf cs format, i.e. :50*ag")
    parser.add_argument("--offset", default=0, help="Amino acid number offset to add to plotting index")
    return parser


def check_required_files(arguments, caller):
    if not os.path.isfile(arguments.library):
        raise Exception("Cannot find library file: " + arguments.library)
    if not os.path.isfile(arguments.reference):
        raise Exception("Cannot find reference file: " + arguments.reference)
    if caller == "classify":
        if not os.path.isfile(arguments.reads):
            raise Exception("Cannot find read file: " + arguments.reads)
    elif caller == "quantify":
        for file in arguments.files:
            if not os.path.isfile(file):
                raise Exception("Cannot find input file: " + file)
    return


def classify():
    arg_parser = build_classify_parser()
    args = arg_parser.parse_args()
    check_required_files(args, "classify")
    analyze(aligned_reads=args.reads, library=args.library, reference=args.reference,
            position=args.position, output=args.output, excluded=args.excluded, filter=args.filter)
    #base_name = "../5G-W2-1M"
    #analyze(aligned_reads=base_name+".paf", library="../dms_libs/5_short.csv", reference="../deaminase.fa",
    #        position=54, output=base_name, excluded="../targets.txt", filter=0)


def quantify():
    arg_parser = build_quantify_parser()
    args = arg_parser.parse_args()
    check_required_files(args, "quantify")
    run_files(sample_name=args.name, filenames=args.files, lib_file=args.library, reference=args.reference,
              pos=args.position, output_dir=args.output, target=args.target, offset=args.offset, save_quants=True)


#if __name__ == "__main__":
#    classify()
