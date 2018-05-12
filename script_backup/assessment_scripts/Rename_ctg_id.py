import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='input file name')

parser.add_argument('-out',
                    required=True,
                    help='input file name')

parser.add_argument('-splitter',
                    required=True,
                    help='splitter')

parser.add_argument('-keep_number',
                    required=True,
                    type=int,
                    help='information want to be kept')

args = vars(parser.parse_args())
input_fasta_file = args['in']
renamed_fasta_file = args['out']
splitter = args['splitter']
keep_number = args['keep_number']


fasta_file_renamed_handle = open(renamed_fasta_file, 'w')
for each_ctg in SeqIO.parse(input_fasta_file, 'fasta'):
    each_ctg_seq = str(each_ctg.seq)
    each_ctg_id = each_ctg.id
    each_ctg_id_split = each_ctg_id.strip().split(splitter)
    each_ctg_id_renamed = ''
    if keep_number == 1:
        each_ctg_id_renamed = '%s' % (each_ctg_id_split[0])
    if keep_number == 2:
        each_ctg_id_renamed = '%s%s%s' % (each_ctg_id_split[0], splitter, each_ctg_id_split[1])
    if keep_number == 3:
        each_ctg_id_renamed = '%s%s%s%s%s' % (each_ctg_id_split[0], splitter, each_ctg_id_split[1], splitter, each_ctg_id_split[2])
    if keep_number == 4:
        each_ctg_id_renamed = '%s%s%s%s%s%s%s' % (each_ctg_id_split[0], splitter, each_ctg_id_split[1], splitter, each_ctg_id_split[2], splitter, each_ctg_id_split[3])
    if keep_number == 5:
        each_ctg_id_renamed = '%s%s%s%s%s%s%s%s%s' % (each_ctg_id_split[0], splitter, each_ctg_id_split[1], splitter, each_ctg_id_split[2], splitter, each_ctg_id_split[3], splitter, each_ctg_id_split[4])

    export_dna_record(each_ctg_seq, each_ctg_id_renamed, '', fasta_file_renamed_handle)

fasta_file_renamed_handle.close()
