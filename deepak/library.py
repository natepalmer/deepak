
class MutationLibrary(dict):
    __slots__ = ["_wt", "file", "reference"]

    def __init__(self, *args, **kwargs):
        super(MutationLibrary, self).__init__(*args, **kwargs)
        self._wt = None
        self.file = None
        self.reference = None

    def construct(self, filename, align_pos):
        """ Initializes empty library from sequence file
        NOTE: generate library treats end mutations as real, while aligners generally soft clip end mutations.
        Ensure amplicon extends at least 30 bp beyond the end of the mutation library region. """
        self.file = filename
        with open(filename) as lib_file:
            for line in lib_file:
                mutation = self.generate_cs(line.strip(), align_pos, normalize=True)
                if mutation != "":
                    self[mutation] = 0
                else:
                    self["wt"] = 0
        return True

    def generate_cs(self, query, alignment_pos, normalize=False):
        """ Generate a paf cs string for a query sequence.
         NOTE: this function cannot handle indels, only substitutions """
        trimmed_ref = self.reference[alignment_pos:]
        cs_fields = list()
        last = 0  # match=0, mismatch=1
        count = 0
        if normalize:
            count -= alignment_pos
        for i, base in enumerate(query):
            try:
                if base != trimmed_ref[i]:
                    if last == 0:
                        cs_fields.append(":{}".format(str(i - count)))
                    cs_fields.append("*{}{}".format(trimmed_ref[i].lower(), base.lower()))
                    last = 1
                    count = i + 1
                else:
                    last = 0
            except IndexError:
                break
        # cs_fields.append(":" + str(i - count))
        return "".join(cs_fields)
