import os
import time
from Bio import AlignIO

def list_to_segments(list_in):

    segments_out = []
    current_element = None
    current_segment = [None, None]
    for each_element in list_in:

        # for the first ellment
        if current_element == None:
            current_element = each_element
            current_segment = [each_element, each_element]

        elif each_element == current_element + 1:
            current_segment[1] = each_element
            current_element = each_element

        elif each_element != current_element + 1:

            # add segment to list
            segments_out.append(current_segment)

            # resetting segment
            current_segment = [each_element, each_element]
            current_element = each_element

    # add segment to list
    segments_out.append(current_segment)


    return segments_out


def slice_string(list_in, col_to_remove):

    all_element_index = list(range(len(list_in)))
    col_to_remove_0_base = [(i - 1) for i in col_to_remove]

    wanted_element_index = []
    for i in all_element_index:
        if i not in col_to_remove_0_base:
            wanted_element_index.append(i)

    # get wanted_segments
    wanted_segments = list_to_segments(wanted_element_index)

    # get concatenated segments
    concatenated_segments_value = []
    for segment in wanted_segments:

        # for single element segment
        if segment[0] == segment[1]:
            segment_value = list_in[segment[0]]

        # for multiple elements segment
        else:
            segment_value = list_in[(segment[0]):(segment[1] + 1)]

        # concatenate segments
        concatenated_segments_value.append(segment_value)

    return concatenated_segments_value


def remove_low_cov_and_consensus_columns(alignment_file_in, minimal_cov, min_consensus, alignment_file_out):

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
            gap_percent = (dash_number / sequence_number) * 100

            if gap_percent > min_cov:
                low_cov_columns.append(n + 1)

            n += 1

        # remove identified columns
        alignment_new = remove_multi_columns_from_msa(alignment, low_cov_columns)
        #alignment_new = slice_string(alignment, low_cov_columns)


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
                each_aa_percent = (current_column.count(each_aa) / sequence_number) * 100
                if each_aa_percent > most_abundant_aa_percent:
                    most_abundant_aa_percent = each_aa_percent

            # if maximum percent lower than provided cutoff, add current column to low consensus column list
            if most_abundant_aa_percent < min_consensus:
                low_css_columns.append(n + 1)

            n += 1

        # remove identified columns
        alignment_new = remove_multi_columns_from_msa(alignment, low_css_columns)
        #alignment_new = slice_string(alignment, low_css_columns)

        return alignment_new


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


wd = '/Users/songweizhi/Desktop/align_format'
os.chdir(wd)

minimal_cov = 50
min_consensus = 25

alignment_file_in = 'TIGR00422_aligned_out.fasta'
alignment_file_out = 'TIGR00422_aligned_out.aln.2'
#alignment_file_in = 'NorthSea_c5_species_tree_tmp.aln'
#alignment_file_out = 'NorthSea_c5_species_tree_tmp.aln.2'


t1 = time.time()
remove_low_cov_and_consensus_columns(alignment_file_in, minimal_cov, min_consensus, alignment_file_out)
t2 = time.time()
print(t2 - t1)

# 6.42

