import os
import re
from time import sleep


def create_append(dictionary, key, value):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]
    return


def resolve_codon(look_ahead_components, location):
    skip = 0
    search_string = ""
    frame = location % 3
    for pos in range(2 - frame):
        if pos == len(look_ahead_components):
            break
        if look_ahead_components[pos][0] in ('-', '+'):
            break
        elif look_ahead_components[pos][0] == ":":
            if (pos == 0) and (int(look_ahead_components[pos][1:]) == 1) \
                    and (len(look_ahead_components) > 1) and (frame == 0):
                if look_ahead_components[1][0] == "*":
                    search_string += ":1"
                    skip += 1
            else:
                break
        elif look_ahead_components[pos][0] == "*":
            search_string += look_ahead_components[pos]
            skip += 1
    return search_string, skip


def chunk_paf(paf_record):
    """ Takes a PafRecord object or a full cs string and returns a list of fields from its cs string """
    if type(paf_record) == str:
        cs = paf_record
    else:
        cs = paf_record.tags["cs"]
    separable_cs = re.sub(r'(?<!\A)([-:*+])', r',\1', cs)
    return separable_cs.split(",")


def make_output_directory(base):
    if base.endswith("/"):
        bn = base[:-1]
    else:
        bn = base.rsplit("/")[-1].split(".")[0]
    out_dir = bn+"/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    else:
        print("Warning: "+out_dir+" already exists, overwriting files in 5 seconds")
        sleep(5)
        print("Overwriting")
    return out_dir
