import os
import shutil
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def get_common_stop_sequence():
    positive_strand_stop_codon = ['TAA', 'TAG', 'TGA']
    negative_strand_stop_codon = ['TTA', 'CTA', 'TCA']
    nc_bases = ['A', 'T', 'G', 'C']
    random_pos = random.sample(positive_strand_stop_codon, 3)
    random_neg = random.sample(negative_strand_stop_codon, 3)
    random_nc = random.sample(nc_bases, 4)
    common_stop_sequence = '%s%s%s%s%s%s%s%s%s%s' % (random_pos[0],
                                                     random_nc[0],
                                                     random_pos[1],
                                                     random_nc[1],
                                                     random_pos[2],
                                                     random_neg[0],
                                                     random_nc[2],
                                                     random_neg[1],
                                                     random_nc[3],
                                                     random_neg[2])
    common_stop_sequence = 'TAGATAAATGATTAGTTAGTTA'
    return common_stop_sequence


def get_random_insertion(recipient_sequence, insert_sequence_list, insert_sequence_id_list):

    insert_gene_number = len(insert_sequence_list)
    length = len(recipient_sequence)
    length_list = list(range(1, length))
    random_bases = random.sample(length_list, insert_gene_number)
    random_bases_sorted = sorted(random_bases)

    # get the start and stop points of all sub_sequences
    sub_sequences_list = []
    n = 0
    first_sequence = [1, random_bases_sorted[n]]
    first_sequence_nc = recipient_sequence[0:random_bases_sorted[n]]
    sub_sequences_list.append(first_sequence)
    while n < insert_gene_number - 1:
        current_sequence = [random_bases_sorted[n] + 1, random_bases_sorted[n+1]]
        sub_sequences_list.append(current_sequence)
        n += 1
    last_sequence = [random_bases_sorted[n] + 1, length]
    sub_sequences_list.append(last_sequence)

    # get new sequences
    new_seq = ''
    m = 0
    while m <= insert_gene_number:
        if m == 0:
            new_seq = first_sequence_nc
        if m > 0:
            current_subsequence_start = sub_sequences_list[m][0] - 1
            current_subsequence_stop = sub_sequences_list[m][1]
            current_subsequence = recipient_sequence[current_subsequence_start:current_subsequence_stop]
            current_insert = insert_sequence_list[m - 1]
            current_insert_id = insert_sequence_id_list[m - 1]
            current_stop_seq_start = get_common_stop_sequence()
            current_stop_seq_end = get_common_stop_sequence()
            current_insert_with_stop = '%s%s%s' % (current_stop_seq_start, current_insert, current_stop_seq_end)
            # print(current_insert_id)
            # print('start: %s' % current_stop_seq_start)
            # print('end: %s' % current_stop_seq_end)
            # print('\n')
            new_seq += current_insert_with_stop
            new_seq += current_subsequence
        m += 1

    return new_seq


parser = argparse.ArgumentParser()

parser.add_argument('-transfer_profile',
                    required=True,
                    help='transfer file')

parser.add_argument('-recipient_genome_folder',
                    required=True,
                    help='recipient genome folder')

parser.add_argument('-combined_ffn_file',
                    required=True,
                    help='combined.ffn file')

parser.add_argument('-recipient_genome_extension',
                    required=False,
                    default='fna',
                    help='extension of recipient genome files')

args = vars(parser.parse_args())

recipients_folder = args['recipient_genome_folder']
if recipients_folder[-1] == '/':
    recipients_folder = recipients_folder[:-1]
combined_ffn_file = args['combined_ffn_file']
transfer_profile_file = args['transfer_profile']

#wd = '/Users/songweizhi/Desktop/GeneTransfer_wd'
wd = os.getcwd()
pwd_combined_ffn = '%s/%s' % (wd, combined_ffn_file)
pwd_transfers = '%s/%s' % (wd, transfer_profile_file)
pwd_recipients_folder = '%s/%s' % (wd, recipients_folder)
recipients_genome_extension = 'fna'
pwd_output_folder = '%s/GeneTransfer_outputs' % wd


if os.path.isdir(pwd_output_folder):
    shutil.rmtree(pwd_output_folder)
    os.mkdir(pwd_output_folder)
else:
    os.mkdir(pwd_output_folder)


# read insert sequences into a list
transfers = open(pwd_transfers)
for each in transfers:
    each_split = each.strip().split(',')
    recipient_genome_id = each_split[0]
    insert_sequence_id_list_unorder = each_split[1:]
    pwd_recipient_genome = '%s/%s.%s' % (pwd_recipients_folder, recipient_genome_id, recipients_genome_extension)
    recipient_genome = SeqIO.read(pwd_recipient_genome, 'fasta')
    recipient_genome_nc = str(recipient_genome.seq)

    # read insert sequences into list
    combined_ffn_handle = SeqIO.parse(pwd_combined_ffn, 'fasta')
    insert_sequence_id_list = []
    insert_sequence_seq_list = []
    for each_gene in combined_ffn_handle:
        if each_gene.id in insert_sequence_id_list_unorder:
            insert_sequence_id_list.append(each_gene.id)
            insert_sequence_seq_list.append(str(each_gene.seq))

    new_seq = get_random_insertion(recipient_genome_nc, insert_sequence_seq_list, insert_sequence_id_list)
    new_seq_record = SeqRecord(Seq(new_seq, IUPAC.unambiguous_dna), id=recipient_genome_id, description='')
    pwd_new_seq_file = '%s/%s.%s' % (pwd_output_folder, recipient_genome_id, recipients_genome_extension)
    new_seq_file_handle = open(pwd_new_seq_file, 'w')
    SeqIO.write(new_seq_record, new_seq_file_handle, 'fasta')
    new_seq_file_handle.close()
