#!/bin/bash

~/Applications/bbmap/reformat.sh in=tests.fa out=tests.fastq overwrite=true
~/Applications/minimap2/minimap2 -x sr --cs reference_test.fa tests.fastq > mapped_tests.paf