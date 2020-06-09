FILTERS = {
    "len_aligned": (400, "ge"),
    "n_mutations": (1, "lt"),
    "indel": (False, "eq"),
    "lib_identity": (1, lambda x, y: len(y) <= x)
}
