import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument('-seq', required=True, help='sequence file')

args = vars(parser.parse_args())
sequence_file = args['seq']

total_length_bp = 0
for each_seq in SeqIO.parse(sequence_file, 'fasta'):
    total_length_bp += len(each_seq.seq)

total_length_Mbp = float("{0:.2f}".format(total_length_bp/(1024*1024)))

print('The total length of sequences in file %s is %s in bp and  %s in Mbp' % (sequence_file, total_length_bp, total_length_Mbp) )

