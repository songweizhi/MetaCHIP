import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def remove_single_columns_from_msa(alignment_in, column_to_remove):

    alignment_column_l = alignment_in[:, :column_to_remove - 1]
    alignment_column_r = alignment_in[:, column_to_remove:]
    alignment_new = alignment_column_l + alignment_column_r

    return alignment_new


def remove_multi_columns_from_msa(alignment_in, column_list):

    alignment_new = alignment_in
    removed_col_num = 0
    for column in sorted(column_list):
        alignment_new = remove_single_columns_from_msa(alignment_new, column - removed_col_num)
        removed_col_num += 1

    return alignment_new


def remove_low_cov_columns(alignment, min_cov):

    # get columns with low coverage
    sequence_number = len(alignment)
    total_col_num = alignment.get_alignment_length()
    low_cov_columns = []
    n = 0
    while n < total_col_num:
        current_column = alignment[:, n]
        dash_number = current_column.count('-')
        gap_percent = (dash_number/sequence_number)*100

        if gap_percent > min_cov:
            low_cov_columns.append(n + 1)

        n += 1

    # remove identified columns
    alignment_new = remove_multi_columns_from_msa(alignment, low_cov_columns)

    return alignment_new


def remove_low_consensus_columns(alignment, min_consensus):

    # get columns with low coverage
    sequence_number = len(alignment)
    total_col_num = alignment.get_alignment_length()
    low_css_columns = []
    n = 0
    while n < total_col_num:
        current_column = alignment[:, n]

        # get all aa in current column
        aa_list = set()
        for aa in current_column:
            aa_list.add(aa)

        # get maximum aa percent
        most_abundant_aa_percent = 0
        for each_aa in aa_list:
            each_aa_percent = (current_column.count(each_aa)/sequence_number)*100
            if each_aa_percent > most_abundant_aa_percent:
                most_abundant_aa_percent = each_aa_percent

        # if maximum percent lower than provided cutoff, add current column to low consensus column list
        if most_abundant_aa_percent < min_consensus:
            low_css_columns.append(n + 1)

        n += 1

    # remove identified columns
    alignment_new = remove_multi_columns_from_msa(alignment, low_css_columns)

    return alignment_new


def remove_low_cov_and_consensus_columns(alignment_file_in, minimal_cov, min_consensus, alignment_file_out):

    # read in alignment
    alignment = AlignIO.read(alignment_file_in, "fasta")

    # remove_low_cov_columns
    alignment_cov = remove_low_cov_columns(alignment, minimal_cov)

    # remove_low_consensus_columns
    alignment_cov_css = remove_low_consensus_columns(alignment_cov, min_consensus)

    # write filtered alignment
    alignment_file_out_handle = open(alignment_file_out, 'w')
    for each_seq in alignment_cov_css:
        alignment_file_out_handle.write('>%s\n' % str(each_seq.id))
        alignment_file_out_handle.write('%s\n' % str(each_seq.seq))
    alignment_file_out_handle.close()


alignment_file_in = '/Users/songweizhi/Desktop/align_format/TIGR00422_aligned_out.fasta'
#alignment_file_in = '/Users/songweizhi/Desktop/align_format/Total2094_g664_species_tree.aln'
minimal_cov = 50
min_consensus = 25
alignment_file_out = '/Users/songweizhi/Desktop/align_format/TIGR00422_aligned_out_filtered.aln'
#remove_low_cov_and_consensus_columns(alignment_file_in, minimal_cov, min_consensus, alignment_file_out)




def remove_multi_columns_from_msa2(alignment_in, column_list):

    alignment_new = alignment_in
    removed_col_num = 0
    for column in sorted(column_list):
        alignment_new = remove_single_columns_from_msa(alignment_new, column - removed_col_num)
        removed_col_num += 1

    return alignment_new


alignment = AlignIO.read(alignment_file_in, "fasta")
column_list = [1, 2, 3, 4, 5, 8, 9, 32, 33, 43, 44, 119, 120, 121, 122, 123, 124, 125, 126, 147, 151, 152, 153, 154, 155, 156, 157, 158]

alignment_new = remove_multi_columns_from_msa2(alignment, column_list)








# write filtered alignment
alignment_file_out_handle = open(alignment_file_out, 'w')
alignment_file_out_handle.write('\n')
alignment_file_out_handle.write('\n')
alignment_file_out_handle.write('\n')

for each_seq in alignment:
    #alignment_file_out_handle.write('>%s\n' % str(each_seq.id))
    alignment_file_out_handle.write('%s\n' % str(each_seq.seq))

alignment_file_out_handle.write('\n')

for each_seq in alignment_new:
    # alignment_file_out_handle.write('>%s\n' % str(each_seq.id))
    alignment_file_out_handle.write('%s\n' % str(each_seq.seq))

alignment_file_out_handle.close()


os.system('cat %s %s > /Users/songweizhi/Desktop/align_format/combined.aln' % (alignment_file_in, alignment_file_out))


######################################## for demo ########################################

# # get the number of sequences in an alignment
# sequence_number = len(alignment)


# # get the number of column
# column_num = alignment_new.get_alignment_length()


# # print the first and last row
# print(alignment[0])
# print(alignment[-1])


# # print the first column
# print(alignment[:, 0])


# # print the first three column
# print(alignment[:, 0:3])
