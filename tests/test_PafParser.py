import sys, os, subprocess
from Bio import SeqIO
import re
import pytest

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from deepak import resolve_codon, PafParser, PafRecord

subprocess.run("./make_tests.sh")
test_object = PafParser("mapped_tests.paf", "reference_test.fa")
test_object.generate_library("library_test.tsv", 27)
print(test_object.library)
fasta_input = list(SeqIO.parse("tests.fa", "fasta"))  # Okay to list for tests as tests.fa will be small
only_subs = [x for x in fasta_input if "+" not in x.name and "-" not in x.name]
look_ahead_tests = {":5": ("", 0), ":1*ag": (":1*ag", 2), ":1+tc": ("", 0), ":1-tc": ("", 0), "*ag:4": ("*ag", 1),
                    "*ag*ag": ("*ag*ag", 2), "*ag+tc": ("*ag", 1), "*ag-tc:4": ("*ag", 1), "+tc:4": ("", 0),
                    "-tc:4": ("", 0)}
scan_fields_answers = {"wild_type": (set(), 0), ":27*ga": (set(), 1), ":27*gc:44*gt*at": ({(":72*gt*at")}, 1),
                       ":72*gt*ag*cg:342*ca*ga*ga*at": ({(":72*gt*ag*cg")}, 2),
                       ":84*ac*gc:78*cg:65*ta:60*ca:123*ct:1*ca*ga*ga*at": (set(), 7),
                       ":72*gt*at:4*ta": ({":72*gt*at", ":78*ta"}, 0),
                       ":27*ga:27*at:361*ct*gt*gt*at": (set(), 4)}
with open("mapped_tests.paf") as paf_in:
    paf_records = [PafRecord(x.strip()) for x in paf_in.readlines()]


def generate_lookahead_answers(question, frame):
    if frame % 3 == 0:
        return look_ahead_tests[question]
    elif frame % 3 == 1 and question[0] == "*":
        return question[:3], 1
    else:
        return "", 0


@pytest.mark.parametrize("sequence_obj", only_subs)
def test_generate_paf(sequence_obj):
    output_cs = test_object.generate_cs(str(sequence_obj.seq), 27, normalize=True)
    if sequence_obj.name == "wild_type":
        assert output_cs == ""
    else:
        assert output_cs == sequence_obj.name


@pytest.mark.parametrize("test", look_ahead_tests.keys())
@pytest.mark.parametrize("frame", range(3, 6))
def test_resolve_codon(test, frame):
    items = re.sub(r'(?<!\A)([-:*+])', r',\1', test).split(",")
    x = resolve_codon(items, frame)
    answer = generate_lookahead_answers(test, frame)
    assert x == answer


@pytest.mark.parametrize("record", paf_records)
def test_scan_fields(record):
    cs = record.tags["cs"]
    result = test_object.scan_fields(record)
    if "-" in cs or "+" in cs:
        assert result[2] == "Indel"
    else:
        assert result[:2] == scan_fields_answers[record.name]
