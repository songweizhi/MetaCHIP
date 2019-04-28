import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
import s0_Kelp_bins_config
import seaborn as sns
import matplotlib.pyplot as plt


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def combine_PG_output(PG_output_file_list_with_path, output_prefix, detection_ranks, combined_PG_output, combined_PG_output_normal):

    HGT_identity_dict = {}
    HGT_end_match_dict = {}
    HGT_full_length_match_dict = {}
    HGT_direction_dict = {}
    HGT_occurence_dict = {}
    HGT_concatenated_list = []
    for pwd_PG_output_file in PG_output_file_list_with_path:
        file_path, file_name = os.path.split(pwd_PG_output_file)
        taxon_rank = file_name[len(output_prefix) + 1]
        if taxon_rank in detection_ranks:
            for PG_HGT in open(pwd_PG_output_file):
                if not PG_HGT.startswith('Gene_1'):
                    PG_HGT_split = PG_HGT.strip().split('\t')

                    gene_1 = PG_HGT_split[0]
                    gene_2 = PG_HGT_split[1]
                    identity = float(PG_HGT_split[4])
                    end_match = PG_HGT_split[5]
                    full_length_match = PG_HGT_split[6]
                    direction = PG_HGT_split[7]
                    concatenated = '%s___%s' % (gene_1, gene_2)

                    if concatenated not in HGT_concatenated_list:
                        HGT_concatenated_list.append(concatenated)

                    # store in dict
                    if concatenated not in HGT_identity_dict:
                        HGT_identity_dict[concatenated] = identity

                    if concatenated not in HGT_end_match_dict:
                        HGT_end_match_dict[concatenated] = end_match

                    if concatenated not in HGT_full_length_match_dict:
                        HGT_full_length_match_dict[concatenated] = full_length_match

                    if direction != 'NA':
                        if concatenated not in HGT_direction_dict:
                            HGT_direction_dict[concatenated] = [direction]
                        else:
                            HGT_direction_dict[concatenated].append(direction)

                    if direction != 'NA':
                        if concatenated not in HGT_occurence_dict:
                            HGT_occurence_dict[concatenated] = [taxon_rank]
                        else:
                            HGT_occurence_dict[concatenated].append(taxon_rank)


    detection_ranks_all = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    detection_ranks_list = []
    for each_rank in detection_ranks_all:
        if each_rank in detection_ranks:
            detection_ranks_list.append(each_rank)


    HGT_occurence_dict_0_1_format = {}
    for each_HGT in HGT_occurence_dict:
        occurence_str = ''
        for each_level in detection_ranks_list:
            if each_level in HGT_occurence_dict[each_HGT]:
                occurence_str += '1'
            else:
                occurence_str += '0'
        HGT_occurence_dict_0_1_format[each_HGT] = occurence_str


    combined_output_handle = open(combined_PG_output, 'w')
    combined_output_handle_normal = open(combined_PG_output_normal, 'w')

    combined_output_handle.write('Gene_1\tGene_2\tIdentity\toccurence(%s)\tend_match\tfull_length_match\tdirection\n' % detection_ranks)
    combined_output_handle_normal.write('Gene_1\tGene_2\tIdentity\toccurence(%s)\tend_match\tfull_length_match\tdirection\n' % detection_ranks)

    for concatenated_HGT in sorted(HGT_concatenated_list):
        concatenated_HGT_split = concatenated_HGT.split('___')

        concatenated_HGT_direction = 'NA'
        if concatenated_HGT in HGT_direction_dict:

            concatenated_HGT_direction_list = HGT_direction_dict[concatenated_HGT]
            concatenated_HGT_direction_list_uniq = unique_list_elements(concatenated_HGT_direction_list)

            if len(concatenated_HGT_direction_list_uniq) == 1:
                concatenated_HGT_direction = concatenated_HGT_direction_list[0]
            else:
                concatenated_HGT_direction = 'both'
                for HGT_direction in concatenated_HGT_direction_list_uniq:
                    HGT_direction_freq = (concatenated_HGT_direction_list.count(HGT_direction))/len(concatenated_HGT_direction_list)
                    if HGT_direction_freq > 0.5:
                        concatenated_HGT_direction = HGT_direction + '(' + str(float("{0:.2f}".format(HGT_direction_freq*100))) + '%)'

        if concatenated_HGT in HGT_occurence_dict_0_1_format:
            occurence_formatted = HGT_occurence_dict_0_1_format[concatenated_HGT]
        else:
            occurence_formatted = '0'*len(detection_ranks_list)

        for_out = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (concatenated_HGT_split[0],
                                                    concatenated_HGT_split[1],
                                                    HGT_identity_dict[concatenated_HGT],
                                                    occurence_formatted,
                                                    HGT_end_match_dict[concatenated_HGT],
                                                    HGT_full_length_match_dict[concatenated_HGT],
                                                    concatenated_HGT_direction)

        combined_output_handle.write(for_out)

        if (HGT_end_match_dict[concatenated_HGT] == 'no') and (HGT_full_length_match_dict[concatenated_HGT] == 'no') and (concatenated_HGT_direction != 'NA') and (concatenated_HGT_direction != 'both'):
            combined_output_handle_normal.write(for_out)


    combined_output_handle.close()
    combined_output_handle_normal.close()


wd = '/Users/songweizhi/Desktop/KelpBins/combined_pcofg/PG_pcofg_new_right_algorithm'
PG_output_folder =          'TT_90MGs_PG'
pwd_candidates_seq_file =   'GoodBins_0.5_0.05_all_combined_ffn.fasta'
output_prefix =             'GoodBins_0.5_0.05'
detection_ranks =           'pcofg'
grouping_file_highest_rank = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_pcofg_grouping/GoodBins_0.5_0.05_p10_grouping.txt'


pwd_candidates_file_PG_txt =            '%s_PG_%s.txt'          % (output_prefix, detection_ranks)
pwd_candidates_file_PG_normal_txt =     '%s_PG_%s_normal.txt'   % (output_prefix, detection_ranks)
pwd_candidates_file_ET_validated_ffn =  '%s_PG_%s_recipient_gene.ffn'   % (output_prefix, detection_ranks)
pwd_candidates_file_ET_validated_faa =  '%s_PG_%s_recipient_gene.faa'   % (output_prefix, detection_ranks)


os.chdir(wd)


PG_output_file_re = '%s/*.txt' % PG_output_folder
PG_output_file_list_with_path = [('%s/%s' % (PG_output_folder, os.path.basename(file_name))) for file_name in glob.glob(PG_output_file_re)]

combine_PG_output(PG_output_file_list_with_path, output_prefix, detection_ranks, pwd_candidates_file_PG_txt, pwd_candidates_file_PG_normal_txt)


#################### export sequence of HGTs as normal ####################

# get recipient gene list
recipient_gene_list = set()
recipient_genome_list = []
for each in open(pwd_candidates_file_PG_normal_txt):
    if not each.startswith('Gene_1'):

        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = each_split[6]

        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        recipient_genome_list.append(recipient_genome)

        if gene_1_genome == recipient_genome:
            recipient_gene_list.add(gene_1)
        if gene_2_genome == recipient_genome:
            recipient_gene_list.add(gene_2)



# export sequence of recipient genes
combined_output_validated_fasta_nc_handle = open(pwd_candidates_file_ET_validated_ffn, 'w')
combined_output_validated_fasta_aa_handle = open(pwd_candidates_file_ET_validated_faa, 'w')
for each_candidate in SeqIO.parse(pwd_candidates_seq_file, 'fasta'):
    if each_candidate.id in recipient_gene_list:
        # output nc sequences
        SeqIO.write(each_candidate, combined_output_validated_fasta_nc_handle, 'fasta')
        # output aa sequences
        each_candidate_aa = each_candidate
        each_candidate_aa.seq = each_candidate_aa.seq.translate()
        SeqIO.write(each_candidate_aa, combined_output_validated_fasta_aa_handle, 'fasta')
combined_output_validated_fasta_nc_handle.close()
combined_output_validated_fasta_aa_handle.close()


########## plot the number of HGT per genome ##########

# get input genome list
input_genome_list = set()
for each_genome in open(grouping_file_highest_rank):
    input_genome_list.add(each_genome.strip().split(',')[1])


# store genome size in dict
genome_size_dict = {}
for each_size in open(s0_Kelp_bins_config.genome_size_file):
    if each_size.strip() != 'Genome\tSize(Mbp)':
        each_size_split = each_size.strip().split('\t')
        genome_file_name = each_size_split[0]
        genome_size_Mbp = float(each_size_split[1])
        genome_size_dict[genome_file_name] = genome_size_Mbp


# get the number of HGT detected from each genome
genome_HGT_num_absolute_dict = {}
for input_genome in input_genome_list:
    genome_HGT_num_absolute_dict[input_genome] = recipient_genome_list.count(input_genome)


genome_HGT_num_norm_dict = {}
for genome in genome_HGT_num_absolute_dict:
    genome_size = genome_size_dict[genome]
    HGT_num_per_Mbp = float("{0:.2f}".format(genome_HGT_num_absolute_dict[genome]/genome_size))
    genome_HGT_num_norm_dict[genome] = HGT_num_per_Mbp



HGT_num_summary_with_label = []
for input_genome in input_genome_list:
    HGT_num_summary_with_label.append([genome_HGT_num_absolute_dict[input_genome], genome_HGT_num_norm_dict[input_genome], genome_size_dict[input_genome], 'HGT'])


HGT_num_summary_with_label_df = pd.DataFrame(np.array(HGT_num_summary_with_label, dtype=object), columns=['absolute', 'norm', 'size', 'label'])
print(HGT_num_summary_with_label_df)



# get plot
sns.set(style='whitegrid')  # whitegrid
ax = sns.boxplot(x='label', y='norm', data=HGT_num_summary_with_label_df, orient='v', width=0.5, showfliers=False)
ax = sns.stripplot(x='label', y='norm', data=HGT_num_summary_with_label_df, color='orange', size=5, orient='v', jitter=0.23)
ax.set(xlabel='Genome', ylabel='Number of HGT/Mbp sequence')
plt.show()





