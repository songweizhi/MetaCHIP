import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='input file')

args = vars(parser.parse_args())
input = args['in']
output_r1 = input.replace('.fasta', '_R1.fasta')
output_r2 = input.replace('.fasta', '_R2.fasta')


r1_handle = open(output_r1, 'w')
r2_handle = open(output_r2, 'w')

n = 0
for each in SeqIO.parse(input, 'fasta'):
    if n % 2 == 0:
        SeqIO.write(each, r1_handle, 'fasta')
    if n % 2 == 1:
        SeqIO.write(each, r2_handle, 'fasta')
    n += 1
r1_handle.close()
r2_handle.close()
