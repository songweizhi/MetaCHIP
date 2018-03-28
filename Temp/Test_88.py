#!/usr/bin/env python
import os
import shutil
import argparse
import itertools
from sys import stdout
from time import sleep
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from string import ascii_uppercase
from Bio.Alphabet import generic_rna

os.chdir('/Users/songweizhi/Desktop/test_translate')

candidates_seq_nc_handle = open('HGT_candidates_nc.fasta', 'w')
candidates_seq_aa_handle = open('HGT_candidates_aa.fasta', 'w')

for each_seq in SeqIO.parse('HGT_candidates.fasta', 'fasta'):

    SeqIO.write(each_seq, candidates_seq_nc_handle, 'fasta')
    each_seq_aa = each_seq
    each_seq_aa.seq = each_seq.seq.translate()
    SeqIO.write(each_seq_aa, candidates_seq_aa_handle, 'fasta')

candidates_seq_nc_handle.close()
candidates_seq_aa_handle.close()




