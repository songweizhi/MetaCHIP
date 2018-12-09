import os
import argparse
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


################################################# input #################################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-HGT_PG', dest='HGT_PG', nargs='?', required=True,  type=str, help='HGT_candidates_PG_validated.txt')
required.add_argument('-grouping', dest='GROUPING', nargs='?', required=True, type=str, help='Grouping file')
required.add_argument('-gt', dest='GT', nargs='?', required=True,  type=str, help='Group to taxon file')
required.add_argument('-size', dest='SIZE', nargs='?', required=True,  type=str, help='Genome size file')
required.add_argument('-prefix', dest='PREFIX', nargs='?', required=True,  type=str, help='Output prefix')

args = vars(parser.parse_args())
HGT_PG_file = args['HGT_PG']
grouping_file = args['GROUPING']
group_to_taxon_file = args['GT']
pwd_genome_size_file = args['SIZE']
output_prefix = args['PREFIX']

################################################# main #################################################

# input files
#wd = '/Users/songweizhi/Desktop/KelpBins/get_HGT_heatmap'
#os.chdir(wd)
# HGT_PG_file =           'HGT_candidates_PG_validated.txt'
# grouping_file =         'GoodBins_0.5_0.05_grouping_c15.txt'
# group_to_taxon_file =   'GoodBins_0.5_0.05_group_to_taxon_c.txt'
# pwd_genome_size_file =  'GoodBins_0.5_0.05_genome_size.txt'
# subset_tree_file = '/Users/songweizhi/Desktop/plot_c/bac120_r86.1.tree'


# output files
HGT_matrix_genome_level_file =              '%s_HGT_matrix_genome.csv'      % output_prefix
HGT_matrix_genome_level_file_normalized =   '%s_HGT_matrix_genome_norm.csv' % output_prefix
HGT_matrix_group_level_file =               '%s_HGT_matrix_group.csv'       % output_prefix
HGT_matrix_group_level_file_normalized =    '%s_HGT_matrix_group_norm.csv'  % output_prefix


############################## store in dict ##############################

# get taxon_name_list for tree subsetting
taxon_name_list = []
group_to_taxon_dict = {}
for each_group in open(group_to_taxon_file):
    taxon_name_list.append(each_group.strip().split(',')[1])
    group_to_taxon_dict[each_group.strip().split(',')[0]] = each_group.strip().split(',')[1]


# store grouping info in dict
genome_to_group_dict = {}
group_to_genome_dict = {}
group_id_list = []
genome_id_list = []
for genome in open(grouping_file):
    genome_split = genome.strip().split(',')

    # get group id list
    if genome_split[0] not in group_id_list:
        group_id_list.append(genome_split[0])

    # get genome id list
    genome_id_list.append(genome_split[1])

    # get genome_to_group_dict
    genome_to_group_dict[genome_split[1]] = genome_split[0]

    # get group_to_genome_dict
    if genome_split[0] not in group_to_genome_dict:
        group_to_genome_dict[genome_split[0]] = [genome_split[1]]
    else:
        group_to_genome_dict[genome_split[0]].append(genome_split[1])


# store genome size in dict
genome_size_dict = {}
for each_size in open(pwd_genome_size_file):
    if each_size.strip() != 'Genome\tSize(Mbp)':
        each_size_split = each_size.strip().split('\t')
        genome_file_name = each_size_split[0]
        genome_file_name_no_ext = '.'.join(genome_file_name.split('.')[:-1])
        genome_size_Mbp = float(each_size_split[1])
        genome_size_dict[genome_file_name_no_ext] = genome_size_Mbp


# get total length of sequence for each group
group_to_length_dict = {}
for each_group in group_to_genome_dict:
    current_genome_member = group_to_genome_dict[each_group]
    group_total_length_Mbp = 0
    for genome in current_genome_member:
        group_total_length_Mbp += genome_size_dict[genome]
    group_to_length_dict[each_group] = group_total_length_Mbp


# read HGT count at genome and group level into dict
HGT_genome_level_dict = {}
HGT_group_level_dict = {}
for each_HGT_pair in open(HGT_PG_file):
    if not each_HGT_pair.startswith('Gene_1'):
        each_HGT_pair_split = each_HGT_pair.strip().split('\t')

        donor_genome = each_HGT_pair_split[6].split('-->')[0]
        donor_group = genome_to_group_dict[donor_genome]
        recipient_genome = each_HGT_pair_split[6].split('-->')[1]
        recipient_group = genome_to_group_dict[recipient_genome]

        donor_to_recipient_genome_level = '%s>%s' % (donor_genome, recipient_genome)
        donor_to_recipient_group_level = '%s>%s' % (donor_group, recipient_group)

        # add to dict genome level
        if donor_to_recipient_genome_level not in HGT_genome_level_dict:
            HGT_genome_level_dict[donor_to_recipient_genome_level] = 1
        else:
            HGT_genome_level_dict[donor_to_recipient_genome_level] += 1

        # add to dict group level
        if donor_to_recipient_group_level not in HGT_group_level_dict:
            HGT_group_level_dict[donor_to_recipient_group_level] = 1
        else:
            HGT_group_level_dict[donor_to_recipient_group_level] += 1


# get matrix (group level)
HGT_matrix_group_level_handle = open(HGT_matrix_group_level_file, 'w')
HGT_matrix_group_level_handle_normalized = open(HGT_matrix_group_level_file_normalized, 'w')
HGT_matrix_group_level_handle.write(',%s\n' % ','.join([group_to_taxon_dict[i] for i in group_id_list]))
HGT_matrix_group_level_handle_normalized.write(',%s\n' % ','.join([group_to_taxon_dict[i] for i in group_id_list]))

for g1 in group_id_list:
    current_donor_value = [group_to_taxon_dict[g1]]
    current_donor_value_norm = [group_to_taxon_dict[g1]]

    for g2 in group_id_list:
        key_g1_g2 = '%s>%s' % (g1, g2)
        current_HGT_num = 0
        if key_g1_g2 in HGT_group_level_dict:
            current_HGT_num = HGT_group_level_dict[key_g1_g2]
        current_donor_value.append(current_HGT_num)
        current_donor_value_norm.append(float("{0:.3f}".format(current_HGT_num/group_to_length_dict[g1])))
    current_donor_value_str = [str(i) for i in current_donor_value]
    current_donor_value_str_norm = [str(i) for i in current_donor_value_norm]
    HGT_matrix_group_level_handle.write('%s\n' % ','.join(current_donor_value_str))
    HGT_matrix_group_level_handle_normalized.write('%s\n' % ','.join(current_donor_value_str_norm))

HGT_matrix_group_level_handle.close()
HGT_matrix_group_level_handle_normalized.close()


# get matrix (genome level)
HGT_matrix_genome_level_handle = open(HGT_matrix_genome_level_file, 'w')
HGT_matrix_genome_level_handle.write(',%s\n' % ','.join(genome_id_list))
for genome_d in genome_id_list:
    current_genome_d_value = [genome_d]
    for genome_r in genome_id_list:
        key_gd_gr = '%s>%s' % (genome_d, genome_r)
        current_gdr_HGT_num = 0
        if key_gd_gr in HGT_genome_level_dict:
            current_gdr_HGT_num = HGT_genome_level_dict[key_gd_gr]
        current_genome_d_value.append(current_gdr_HGT_num)
    current_genome_d_value_str = [str(i) for i in current_genome_d_value]
    HGT_matrix_genome_level_handle.write('%s\n' % ','.join(current_genome_d_value_str))
HGT_matrix_genome_level_handle.close()

