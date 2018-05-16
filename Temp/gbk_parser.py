#!/usr/bin/env python
import sys
import os
from Bio import SeqIO

__author__ = 'weizhisong'

usage = """
####################################################################

  Usage:
  python gbk_parser.py aaa.gbk keep_list.txt

  Keep list file format: One gene name per line.

####################################################################
"""

if len(sys.argv) != 3:
    print(usage)
    exit(1)


# Generate keep list:
keep_list = sys.argv[2]
matches = open(keep_list)
match_list = []
for match in matches:
    match = match.strip()
    match_list.append(match)


# Parse gbk file
input_gbk_file = sys.argv[1]
os.system('mkdir output_gbk')
records = SeqIO.parse(input_gbk_file, 'genbank')
for record in records:
    for feature in record.features:
        if 'locus_tag' in feature.qualifiers:
            for each_mach in match_list:
                if each_mach in feature.qualifiers["locus_tag"]:
                    print(record)
                    SeqIO.write(record, './output_gbk/' + each_mach + '.gbk', 'genbank')
                    SeqIO.write(record, './output_gbk/' + each_mach + '.fasta', 'fasta')
