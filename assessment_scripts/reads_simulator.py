import os
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_header()
    fasta_out.write_record(seq_record)
    fasta_out.write_footer()


def simulate_reads(pwd_genome_file, read_number, read_length, insert_size, output_folder):
    path, file_name = os.path.split(pwd_genome_file)
    genome_name, ext = os.path.splitext(file_name)

    if args['split'] == True:
        output_r1 = '%s/%s_R1.fasta' % (output_folder, genome_name)
        output_r2 = '%s/%s_R2.fasta' % (output_folder, genome_name)
        output_r1_handle = open(output_r1, 'w')
        output_r2_handle = open(output_r2, 'w')
    else:
        output_combined = '%s/%s_R1_R2.fasta' % (output_folder, genome_name)
        output_combined_handle = open(output_combined, 'w')

    seq = str(SeqIO.read(pwd_genome_file, 'fasta').seq)
    sequence_length = len(seq)
    fragment_length = 2 * read_length + insert_size

    n = 1
    while n <= read_number:
        rdm_num = random.randint(1, sequence_length)
        current_fragment = ''
        if (rdm_num + fragment_length) <= sequence_length:
            current_fragment = seq[rdm_num - 1: rdm_num + fragment_length - 1]
        elif (rdm_num + fragment_length) > sequence_length:
            seq_part_1_seq = seq[rdm_num - 1:]
            seq_part_2_seq = seq[:fragment_length - sequence_length + rdm_num - 1]
            current_fragment = seq_part_1_seq + seq_part_2_seq
        current_fragment_r1 = current_fragment[:read_length]
        current_fragment_r2 = current_fragment[-read_length:]
        current_fragment_r2_reverse_complement = str(Seq(current_fragment_r2, generic_dna).reverse_complement())
        current_read_r1_id = 'r%s_from_%s_%sth_bp_#0/1' % (n, genome_name, rdm_num)
        current_read_r2_id = 'r%s_from_%s_%sth_bp_#0/2' % (n, genome_name, rdm_num)
        if args['split'] == 1:
            export_dna_record(current_fragment_r1, current_read_r1_id, '', output_r1_handle)
            export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_r2_handle)
        else:
            export_dna_record(current_fragment_r1, current_read_r1_id, '', output_combined_handle)
            export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_combined_handle)
        # report process
        if (n/5000).is_integer() == 1:
            print('Simulated %s reads' % n)
        n += 1
    if args['split'] == 1:
        output_r1_handle.close()
        output_r2_handle.close()
    else:
        output_combined_handle.close()


####################################################### arguments ######################################################

parser = argparse.ArgumentParser()

parser.add_argument('-r',
                    required=True,
                    help='reference genome')

parser.add_argument('-o',
                    required=True,
                    help='output folder')

parser.add_argument('-n',
                    required=True,
                    type=int,
                    help='reads number')

parser.add_argument('-l',
                    required=True,
                    type=int,
                    help='reads length')

parser.add_argument('-i',
                    required=True,
                    type=int,
                    help='insert size')

parser.add_argument('-split',
                    action="store_true",
                    help='Export forward and reverse reads to separate files')

args = vars(parser.parse_args())
pwd_genome_file = args['r']
output_folder = args['o']
read_number = args['n']
read_length = args['l']
insert_size = args['i']

####################################################### main code ######################################################

simulate_reads(pwd_genome_file, read_number, read_length, insert_size, output_folder)
