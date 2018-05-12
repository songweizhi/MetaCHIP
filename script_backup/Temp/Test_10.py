import os
import shutil
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def get_genome_size(fasta_file):
    genome = SeqIO.parse(fasta_file, 'fasta')
    total_length = 0
    contig_num = 0
    contig_list = []
    contig_length_list = []
    contig_seq_list = []
    for each_contig in genome:
        total_length += len(each_contig.seq)
        contig_num += 1
        contig_list.append(each_contig.id)
        contig_length_list.append(len(each_contig.seq))
        contig_seq_list.append(each_contig.seq)
    return total_length, contig_num, contig_list, contig_length_list, contig_seq_list


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_header()
    fasta_out.write_record(seq_record)
    fasta_out.write_footer()


####################################################### arguments ######################################################

parser = argparse.ArgumentParser()

parser.add_argument('-r',
                    required=True,
                    help='reference genomes')

parser.add_argument('-n',
                    required=True,
                    type=int,
                    help='reads number')

parser.add_argument('-l',
                    required=True,
                    type=int,
                    help='reads length')

parser.add_argument('-c',
                    required=True,
                    type=int,
                    help='circulized genome (1) or not (0)')

parser.add_argument('-o',
                    required=True,
                    help='output file')

args = vars(parser.parse_args())
genome_file = args['r']
reads_number = args['n']
reads_length = args['l']
circulized = args['c']
output_file = args['o']

####################################################### main code ######################################################


wd = os.getcwd()
genome_name, genome_ext = os.path.splitext(genome_file)

# get basics statistics of input genome
genome_size, contig_number, contig_list, contig_length_list, contig_seq_list= get_genome_size(genome_file)

#print('genome_size: %sbp' % genome_size)
#print('contig_number: %s' % contig_number)
#print('contig_list')
#print(contig_list)
#print(contig_seq_list)

reads_number_list = []
weight_per_base = reads_number/genome_size
for each_length in contig_length_list:
    current_reads_number = each_length * weight_per_base
    current_reads_number = int(round(current_reads_number, 0))
    reads_number_list.append(current_reads_number)

#print('reads_number_list')
#print(reads_number_list)

output_file_handle = open(output_file, 'w')
n = 0
ctg_num = 1
while n < contig_number:
    current_ctg_id = contig_list[n]
    current_ctg_seq = contig_seq_list[n]
    current_ctg_seq_length = len(current_ctg_seq)
    current_ctg_reads_allocation = reads_number_list[n]
    #print(current_ctg_id)
    #print(current_ctg_reads_allocation)
    #print(current_ctg_seq)

    m = 1
    while m <= current_ctg_reads_allocation:
        rdm_num = random.randint(1, current_ctg_seq_length)
        current_read_seq = ''
        if rdm_num + reads_length <= current_ctg_seq_length:
            current_read_seq = current_ctg_seq[rdm_num - 1:rdm_num + reads_length - 1]
        if rdm_num + reads_length > current_ctg_seq_length:
            if circulized == 1:
                seq_part_1_seq = current_ctg_seq[rdm_num - 1:]
                seq_part_2_seq = current_ctg_seq[: reads_length - current_ctg_seq_length + rdm_num - 1]
                current_read_seq = seq_part_1_seq + seq_part_2_seq
            if circulized == 0:
                current_read_seq = current_ctg_seq[rdm_num - 1:]

        current_read_id = '%s_%s_r%s' % (genome_name, current_ctg_id, ctg_num)
        if len(current_read_seq) >= 50:
            export_dna_record(str(current_read_seq), current_read_id, '', output_file_handle)
        m += 1
        ctg_num += 1
    n += 1
output_file_handle.close()
