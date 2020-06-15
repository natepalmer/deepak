from deepak.utilities import chunk_paf


class PafRecord:
    __slots__ = ["name", "length", "start", "end", "strand", "ref_name", "ref_length", "ref_start", "ref_end",
                 "len_matches", "len_aligned", "mapping_qual", "tags", "n_mutations", "lib_identity", "indel"]

    def __init__(self, text):
        fields = text.split(sep="\t")
        if len(fields) < 12:
            raise Exception("Input string is missing required components; see pairwise mapping format specification")
        self.name = fields[0]
        self.length = int(fields[1])
        self.start = int(fields[2])
        self.end = int(fields[3])
        self.strand = fields[4]
        self.ref_name = fields[5]
        self.ref_length = int(fields[6])
        self.ref_start = int(fields[7])
        self.ref_end = int(fields[8])
        self.len_matches = int(fields[9])
        self.len_aligned = int(fields[10])
        self.mapping_qual = int(fields[11])
        self.tags = dict()
        for x in fields[12:]:
            y = x.partition(":")
            self.tags[y[0]] = y[2].partition(":")[2]
        #self.tags = {x.partition(":")[0]: x.partition(":")[2].partition(":")[2] for x in fields[12:]}
        self.n_mutations = None
        self.lib_identity = None
        self.indel = False
        self.normalize_position()

    def add_alignment_info(self, lib_identity=None, n_mutations=None, replace=False):
        """
        Add or update information about library membership and/or mutation number.

        *lib_identity* is a set containing the library member identifier strings
        *n_mutations* is an integer specifying the number of additional mutations to add to the count
        """
        if n_mutations is not None:
            if self.n_mutations is None or replace:
                self.n_mutations = n_mutations
            else:
                n_mutations += n_mutations
        if lib_identity is not None:
            if self.lib_identity is None or replace:
                self.lib_identity = lib_identity
            else:
                self.lib_identity.update(lib_identity)

    def normalize_position(self):
        if self.ref_start > 0:
            cs_chunks = chunk_paf(self)
            if len(cs_chunks) == 0:
                return
            elif not cs_chunks[0].startswith(":"):
                self.tags["cs"] = ":{}{}".format(str(self.ref_start), self.tags["cs"])
            else:
                cs_chunks[0] = ":{}".format(str(self.ref_start+int(cs_chunks[0][1:])))
                self.tags["cs"] = "".join(cs_chunks)
        return

    def to_csv(self):
        return [self.name, self.length, self.start, self.end, self.ref_start, self.ref_end, self.len_matches,
                self.len_aligned, self.tags["cs"], self.n_mutations, self.str_lib()]

    def str_lib(self):
        if self.lib_identity is None:
            return "None"
        else:
            return "|".join(self.lib_identity)

    def __len__(self):
        return self.length

    def __str__(self):
        return "PAF record {0} aligned to {1} at pos {2} with quality {3}".format(
            self.name, self.ref_name, str(self.ref_start), str(self.mapping_qual))

    def __repr__(self):
        return "_".join([self.name, self.tags["cs"]])
