COLS = [
    "name", "length", "start", "end", "ref_start", "ref_end", "len_matches", "len_aligned",
    "cs_tag", "n_mutations", "lib_identity"
]

#: List of amino acids in row order for sequence-function maps.
AA_LIST = [
    'H', 'K', 'R',                 # (+)
    'D', 'E',                      # (-)
    'C', 'M', 'N', 'Q', 'S', 'T',  # Polar-neutral
    'A', 'I', 'L', 'V',            # Non-polar
    'F', 'W', 'Y',                 # Aromatic
    'G', 'P'                       # Unique
]
#   '*']
