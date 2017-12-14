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
    for each_contig in genome:
        total_length += len(each_contig.seq)
    return total_length


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_header()
    fasta_out.write_record(seq_record)
    fasta_out.write_footer()


def simulate_reads(pwd_genome_file, read_number, read_length, insert_size, split, output_folder):
    path, file_name = os.path.split(pwd_genome_file)
    genome_name, ext = os.path.splitext(file_name)

    if split == 1:
        output_r1 = '%s/%s_R1.fasta' % (output_folder, genome_name)
        output_r2 = '%s/%s_R2.fasta' % (output_folder, genome_name)
        output_r1_handle = open(output_r1, 'w')
        output_r2_handle = open(output_r2, 'w')
    else:
        output_combined = '%s/%s_R12.fasta' % (output_folder, genome_name)
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
        if split == 1:
            export_dna_record(current_fragment_r1, current_read_r1_id, '', output_r1_handle)
            export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_r2_handle)
        else:
            export_dna_record(current_fragment_r1, current_read_r1_id, '', output_combined_handle)
            export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_combined_handle)
        n += 1
    if split == 1:
        output_r1_handle.close()
        output_r2_handle.close()
    else:
        output_combined_handle.close()


####################################################### arguments ######################################################


parser = argparse.ArgumentParser()

parser.add_argument('-R',
                    required=True,
                    help='reference genomes')

parser.add_argument('-a',
                    required=True,
                    help='abundance file')

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

parser.add_argument('-o',
                    required=True,
                    help='output folder')

args = vars(parser.parse_args())
genome_folder = args['R']
abundance_file = args['a']
reads_number = args['n']
reads_length = args['l']
insert_size = args['i']
split = args['split']
output_folder = args['o']
SimSeq_wd = output_folder

####################################################### main code ######################################################

# prepare folders
wd = os.getcwd()
pwd_genome_folder = '%s/%s' % (wd, genome_folder)
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_qsub_file_folder = '%s/qsub_files' % wd
pwd_SimSeq_wd = '%s/%s' % (wd, SimSeq_wd)

# get genome size list
abundance_combined = open(pwd_abundance_file)
genome_name_list = []
genome_name_no_extension_list = []
genome_size_list = []
replicate_number = 0
for each in open(pwd_abundance_file):
    each_split = each.strip().split('\t')
    genome_name = each_split[0]
    genome_name_no_extension = genome_name.split('.')[0]
    genome_file = '%s/%s/%s' % (wd, genome_folder, genome_name)
    genome_size = get_genome_size(genome_file)
    genome_size_mbp = float("{0:.2f}".format(genome_size / (1024 * 1024)))
    replicate_number = len(each_split) - 1
    genome_name_no_extension_list.append(genome_name_no_extension)
    genome_name_list.append(genome_name)
    genome_size_list.append(genome_size_mbp)

# prepare GemSIM working directory
# if os.path.isdir(pwd_SimSeq_wd):
#     shutil.rmtree(pwd_SimSeq_wd)
#     os.mkdir(pwd_SimSeq_wd)
# else:
#     os.mkdir(pwd_SimSeq_wd)

# if os.path.isdir(pwd_qsub_file_folder):
#     shutil.rmtree(pwd_qsub_file_folder)
#     os.mkdir(pwd_qsub_file_folder)
# else:
#     os.mkdir(pwd_qsub_file_folder)

# main
n = 1
while n <= replicate_number:
    current_abundance_list = []
    for each in open(pwd_abundance_file):
        each_split = each.strip().split('\t')
        current_abundance = int(each_split[n])
        current_abundance_list.append(current_abundance)

    m = 0
    weight_list = []
    while m < len(genome_name_no_extension_list):
        current_weight = (genome_size_list[m])*(current_abundance_list[m])
        weight_list.append(current_weight)
        m += 1
    weight_summation = sum(weight_list)
    multiply_unit = reads_number/weight_summation
    reads_number_list = []
    for each in weight_list:
        current_reads_number = multiply_unit * each
        reads_number_list.append(int(current_reads_number))

    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_SimSeq_wd, n)
    os.mkdir(pwd_current_replicate_wd)

    p = 0
    for current_genome in genome_name_list:
        current_genome_no_extension = current_genome.split('.')[0]
        pwd_current_genome = '%s/%s/%s' % (wd, genome_folder, current_genome)
        print('Processing reapicate_%s, simulating %s reads from %s' % (n, reads_number_list[p], current_genome))
        simulate_reads(pwd_current_genome, reads_number_list[p], reads_length, insert_size, split, pwd_current_replicate_wd)
        # current_wd = os.getcwd()
        # os.chdir(pwd_qsub_file_folder)
        # os.chdir(current_wd)
        p += 1

    # combined all reads together


    if split == 1:
        os.system('cat %s/*_R1.fasta > %s/combined_R1.fasta' % (pwd_current_replicate_wd, output_folder))
        os.system('cat %s/*_R2.fasta > %s/combined_R2.fasta' % (pwd_current_replicate_wd, output_folder))
    else:
        os.system('cat %s/*_R12.fasta > %s/combined.fasta' % (pwd_current_replicate_wd, output_folder))

    os.system('rm -r %s' % pwd_current_replicate_wd)

    n += 1


