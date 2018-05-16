"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
import sys
import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from BCBio import GFF

def main(gff_file, fasta_file):
    out_file = "%s.gbk" % os.path.splitext(gff_file)[0]
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(gff_iter, out_file, "genbank")


gff_file = '/Users/songweizhi/Desktop/000/NorthSea/NorthSea_bin014.gff'
fa_file = '/Users/songweizhi/Desktop/000/NorthSea/NorthSea_bin014.fasta'

os.chdir('/Users/songweizhi/Desktop/000/NorthSea')
main(gff_file, fa_file)




