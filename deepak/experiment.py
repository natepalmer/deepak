import csv
import os
from collections import Counter
from io import StringIO

import pandas as pd
from Bio import SeqIO

from deepak.utilities import chunk_paf, resolve_codon
from deepak.pafparser import PafRecord
from deepak.globals import COLS
from deepak.configuration import FILTERS
from deepak.library import MutationLibrary
from deepak.refiner import Filters


class Experiment:

    def __init__(self, file_name, reference):
        self.n_records = 0
        self.categories = ("Alignment", "Indels", "Many mutations", "Multiple identity", "Valid")
        self.data_files = dict()
        self.filtered_reads = dict()
        self.read_numbers = {key: 0 for key in self.categories}
        self.misaligned = list()
        self.fn = file_name
        self.reference = str(SeqIO.read(reference, "fasta").seq)
        self.library = MutationLibrary()
        self.library.reference = self.reference
        self.target_mutations = list()
        self.filters = Filters()
        self.mutations = Counter()
        self.set_filters()

    def create_output_structure(self, out_dir):
        read_data = "{}data".format(out_dir)
        if not os.path.isdir(read_data):
            os.mkdir(read_data)
        for cat in self.categories:
            file = open(f"{read_data}/{cat}.csv", mode="w")
            writer = csv.writer(file)
            writer.writerow(COLS)  # header
            self.filtered_reads[cat] = (writer, file)
        return

    def initialize_writers(self):
        return {key: self.make_writer() for key in self.categories}

    def make_writer(self):
        stream = StringIO()
        writer = csv.writer(stream)
        return writer, stream

    def set_filters(self, tests=FILTERS):
        for k, v in tests.items():
            assert type(k) is str
            self.filters[k] = v
        return

    def set_excluded_mutations(self, mutation_list):
        self.target_mutations = mutation_list

    def note_failed_alignment(self, query):
        #self.filtered_reads["Alignment"][0].writerow(query.to_csv())
        #if query.name not in self.misaligned:
        #self.misaligned.append(query.name)
        self.read_numbers["Alignment"] += 1
        return

    def common_mutations(self):
        return self.mutations.most_common()

    def scan_fields(self, x):
        components = chunk_paf(x)
        location = 0
        skip = 0
        n_mutations = 0
        filter_failed = None
        lib_identity = set()
        for i, field in enumerate(components):
            if skip > 0:
                skip -= 1
                continue
            if field[0] == ":":
                location += int(field[1:])
            elif field[0] == "*":
                search_string = field
                if i != 0:
                    search_string = ":"+str(location)+search_string
                if i != len(components) - 1:
                    add, skip = resolve_codon(components[i + 1:], location)
                    search_string += add
                    location += skip
                if search_string in self.library:
                    lib_identity.add(search_string)
                else:
                    if search_string not in self.target_mutations:
                        n_mutations += 1
                    self.mutations[search_string] += 1
                location += 1
            else:  # Indel
                filter_failed = "Indel"
                n_mutations += 1
                self.mutations[":" + str(location) + field] += 1
                if field[0] == "-":
                    location += len(field[1:])
        return lib_identity, n_mutations, filter_failed

    def set_read_status(self, paf_record, category):
        if category == "Valid":
            res = paf_record.to_csv()
            self.filtered_reads[category][0].writerow(res)
        self.read_numbers[category] += 1
        return

    def classify(self, report_number=100000):
        with open(self.fn) as paf_file:
            for i, record in enumerate(paf_file):
                if i % report_number == 0:
                    print("Processed {} reads, {} passed filtering".format(i, self.read_numbers["Valid"]), flush=True)
                self.n_records += 1
                x = PafRecord(record.strip())
                if self.filters.test(x, "len_aligned")[0]:
                    lib_identity, n_mutations, filter_failed = self.scan_fields(x)
                    x.add_alignment_info(lib_identity, n_mutations)
                    matches = len(lib_identity)
                    if filter_failed == "Indel" and self.filters["indel"]:
                        self.set_read_status(x, "Indels")
                    elif n_mutations > self.filters["substitution"]:
                        self.set_read_status(x, "Many mutations")
                    elif matches > 1:
                        self.set_read_status(x, "Multiple identity")
                    elif matches in (0, 1):
                        if matches == 0:
                            self.library["wt"] += 1
                            lib_identity.add("wt")
                        self.set_read_status(x, "Valid")
                        if matches == 1:
                            self.library[lib_identity.pop()] += 1
                else:
                    self.note_failed_alignment(x)
        # close csv files
        for cat, double in self.filtered_reads.items():
            double[1].close()
        return True

    def export_summary(self, outer_dir):
        summary = "{}summary".format(outer_dir)
        if not os.path.isdir(summary):
            os.mkdir(summary)
        with open(summary + "/lib.tsv", mode='w') as lib_out:
            print("\n".join(['{}\t{}'.format(k, v) for k, v in self.library.items()]), file=lib_out)
        with open(summary + "/common_mutations.tsv", mode="w") as mut_out:
            print("\n".join(["{}\t{}".format(k[0], k[1]) for k in self.common_mutations()]), file=mut_out)
        with open(summary + "/stats.tsv", mode="w") as stats_out:
            print("\n".join(["{}\t{}".format(k, self.read_numbers[k]) for k in self.read_numbers]), file=stats_out)
        return

    def run(self, out_dir):
        self.create_output_structure(out_dir)
        self.classify()
        self.export_summary(out_dir)
        return

    def export_data(self, outer_dir):
        read_data = "{}data".format(outer_dir)
        if not os.path.isdir(read_data):
            os.mkdir(read_data)
        for cat in self.categories:
            self.filtered_reads[cat].to_csv("{}/{}.csv".format(read_data, cat))
        return

    def export(self, outer_dir):
        self.export_summary(outer_dir)
        self.export_data(outer_dir)
        return
