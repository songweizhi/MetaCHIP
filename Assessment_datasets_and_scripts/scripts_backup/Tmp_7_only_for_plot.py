import os
import re
import sys
import glob
import shutil
import argparse
import numpy as np
from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def plot_identity_distribution(identity_list, plot_file):
    num_bins = 50
    plt.hist(identity_list,
             num_bins,
             alpha=0.1,
             normed=0,  # normed = 1 normalized to 1, that is probablity
             facecolor='blue')
    plt.title('Identity distribution of identified %s HGTs' % len(identity_list))
    plt.xlabel('Identity (%)')
    plt.ylabel('Number of identified HGT')
    plt.subplots_adjust(left=0.15)

    # Get plot
    plt.savefig(plot_file, dpi=300)
    plt.close()


# input
wd = '/Users/songweizhi/Desktop/test_plot'
os.chdir(wd)

pwd_candidates_file_ET =            'in/HGT_candidates_PG.txt'
pwd_candidates_file_ET_validated =  'in/HGT_candidates_PG_validated.txt'
pwd_grouping_file =                 'in/GoodBins_0.5_0.05_grouping_o34.txt'
pwd_genome_size_file =              'in/GoodBins_0.5_0.05_genome_size.txt'
pwd_grouping_id_to_taxon_file =     'in/GoodBins_0.5_0.05_group_to_taxon_o.txt'
circos_HGT_R =                      '/Users/songweizhi/PycharmProjects/MetaCHIP/circos_HGT.R'

# output
pwd_plot_at_ends_number =                           'GoodBins_0.5_0.05_plot_at_ends_stat.png'
pwd_plot_identity_distribution_BM =                 'GoodBins_0.5_0.05_plot_HGT_identities_BM.png'
pwd_plot_identity_distribution_PG =                 'GoodBins_0.5_0.05_plot_HGT_identities_PG.png'
pwd_candidates_file_ET_validated_STAT_genome_txt =  'HGT_candidates_PG_validated_genome_stats.txt'
pwd_candidates_file_ET_validated_STAT_group_txt =   'HGT_candidates_PG_validated_group_stats.txt'
pwd_candidates_file_ET_validated_STAT_png =         'HGT_candidates_PG_validated_stats.png'
pwd_plot_circos =                                   'GoodBins_0.5_0.05_plot_circos_PG.png'


########################################################################################################################
######################################################## Get_plot ######################################################
########################################################################################################################

# read in prediction results
HGT_num_BM_at_end = 0
HGT_num_PG_at_end = 0
HGT_num_BM_not_at_end = 0
HGT_num_PG_not_at_end = 0
identity_list_not_end_BM = []
identity_list_not_end_PG = []
for each_HGT in open(pwd_candidates_file_ET):
    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')
        identity = float(each_HGT_split[4])
        at_end = each_HGT_split[5]
        PG_validation = each_HGT_split[6]
        if at_end == 'no':
            identity_list_not_end_BM.append(identity)
            HGT_num_BM_not_at_end += 1

            if PG_validation != 'N/A':
                identity_list_not_end_PG.append(identity)
                HGT_num_PG_not_at_end += 1

        else: # at_end == 'yes':
            HGT_num_BM_at_end += 1
            if PG_validation != 'N/A':
                HGT_num_PG_at_end += 1


################################################## plot at_ends_stat ###################################################

n_groups = 2
not_at_end_list = (HGT_num_BM_not_at_end, HGT_num_PG_not_at_end)
at_end_list = (HGT_num_BM_at_end, HGT_num_PG_at_end)

# create plot
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.25

rects1 = plt.bar(index, not_at_end_list, bar_width, alpha=0.4, color='g', label='Not at ends', align='center')
rects2 = plt.bar(index + bar_width, at_end_list, bar_width, alpha=0.4, color='b', label='At ends', align='center')

plt.ylabel('Number of predicted HGTs')
plt.title('Location of predicted HGTs')
plt.xticks(index + bar_width/2, ('Best-match', 'Phylogenetic'))
plt.legend()

plt.tight_layout()
plt.savefig(pwd_plot_at_ends_number, dpi=300)
plt.close()


#################################### plot identity distribution of identified HGTs #####################################

# for predicted HGTs:
# 1. not at ends

# plot Identity distribution of identified HGTs
plot_identity_distribution(identity_list_not_end_BM, pwd_plot_identity_distribution_BM)
plot_identity_distribution(identity_list_not_end_PG, pwd_plot_identity_distribution_PG)


################################ plot number of HGT detected from each genome and group ################################

# for predicted HGTs:
# 1. not at ends
# 2. PG validated


# get genome member in each group
group_to_genome_dict = {}
for each_grouping in open(pwd_grouping_file):
    each_grouping_split = each_grouping.strip().split(',')
    group_id = each_grouping_split[0]
    genome_id = each_grouping_split[1]
    if group_id not in group_to_genome_dict:
        group_to_genome_dict[group_id] = [genome_id]
    else:
        group_to_genome_dict[group_id].append(genome_id)


# store genome size in dict
genome_size_dict = {}
for each_size in open(pwd_genome_size_file):
    if each_size.strip() != 'Genome\tSize(Mbp)':
        each_size_split = each_size.strip().split('\t')
        genome_file_name = each_size_split[0]
        genome_file_name_no_ext = '.'.join(genome_file_name.split('.')[:-1])
        genome_size_Mbp = float(each_size_split[1])
        genome_size_dict[genome_file_name_no_ext] = genome_size_Mbp


# get genome_to_HgtNum_dict
genome_to_HgtNum_dict = {}
for each_pair in open(pwd_candidates_file_ET_validated):
    if not each_pair.startswith('Gene_1'):

        each_pair_split = each_pair.strip().split('\t')
        gene_1 = each_pair_split[0]
        gene_2 = each_pair_split[1]
        genome_1 = '_'.join(gene_1.split('_')[:-1])
        genome_2 = '_'.join(gene_2.split('_')[:-1])

        if genome_1 not in genome_to_HgtNum_dict:
            genome_to_HgtNum_dict[genome_1] = 1
        else:
            genome_to_HgtNum_dict[genome_1] += 1

        if genome_2 not in genome_to_HgtNum_dict:
            genome_to_HgtNum_dict[genome_2] = 1
        else:
            genome_to_HgtNum_dict[genome_2] += 1


# get the number of HGT in each group
group_list = []
for each_HGT_pair in open(pwd_candidates_file_ET_validated):
    if not each_HGT_pair.startswith('Gene_1'):
        each_HGT_pair_split = each_HGT_pair.strip().split()
        group_list += [each_HGT_pair_split[2], each_HGT_pair_split[3]]
group_list_uniq = sorted(unique_list_elements(group_list))
group_list_uniq_count = [group_list.count(i) for i in group_list_uniq]


# get total length of sequence for each group
group_to_length_dict = {}
for each_group in group_to_genome_dict:
    current_genome_member = group_to_genome_dict[each_group]
    group_total_length_Mbp = 0
    for genome in current_genome_member:
        group_total_length_Mbp += genome_size_dict[genome]
    group_to_length_dict[each_group] = group_total_length_Mbp


# read group_2_taxon into dict
group_2_taxon_dict = {}
for each_group_2_taxon in open(pwd_grouping_id_to_taxon_file):
    each_group_2_taxon_split = each_group_2_taxon.strip().split(',')
    group_2_taxon_dict[each_group_2_taxon_split[0]] = each_group_2_taxon_split[1]


group_id_with_taxon = []
for each_group in group_list_uniq:
    each_group_new = '(%s) %s' % (each_group, group_2_taxon_dict[each_group])
    group_id_with_taxon.append(each_group_new)


# get group_id_list
group_id_list = []
for each_group in group_to_genome_dict:
    group_id_list.append(each_group)


# normalize HGT number per group with total length
group_list_uniq_count_normalized = [float("{0:.2f}".format(group_list.count(i)/group_to_length_dict[i])) for i in group_list_uniq]


# define color list
color_list = ['black', 'silver', 'darksalmon', 'blueviolet', 'sienna', 'tan', 'purple', 'gold', 'palegreen', 'paleturquoise', 'slategray', 'royalblue', 'plum', 'olivedrab', 'seagreen', 'darkorchid', 'darkkhaki']*100

output_txt_handle = open(pwd_candidates_file_ET_validated_STAT_genome_txt, 'w')
output_txt_handle.write('Group\tGenome\tSize\tHGT\tHGT/Mbp\n')
genome_list_according_group = []
genome_list_according_group_with_group = []
num_list_according_group = []
num_list_according_group_norm = []
color_list_according_group = []
color_index = 1
for each_group_id in sorted(group_id_list):
    current_group_genomes = group_to_genome_dict[each_group_id]
    for each_genome in current_group_genomes:
        each_genome_with_group = '(%s)%s' % (each_group_id, each_genome)
        if each_genome in genome_to_HgtNum_dict:
            genome_list_according_group.append(each_genome)
            genome_list_according_group_with_group.append(each_genome_with_group)

            num_HGT = genome_to_HgtNum_dict[each_genome]
            num_HGT_norm = num_HGT/genome_size_dict[each_genome]
            num_HGT_norm = float("{0:.2f}".format(num_HGT_norm))

            num_list_according_group.append(num_HGT)
            num_list_according_group_norm.append(num_HGT_norm)

            color_list_according_group.append(color_list[color_index])
            for_out = '%s\t%s\t%s\t%s\t%s\n' % (each_group_id, each_genome, genome_size_dict[each_genome], num_HGT, num_HGT_norm)
            output_txt_handle.write(for_out)
    color_index += 1
output_txt_handle.close()


# write stats to file
HGT_PG_STAT_handle = open(pwd_candidates_file_ET_validated_STAT_group_txt, 'w')
n = 0
HGT_PG_STAT_handle.write('Group\tSize(Mbp)\tHGT\tHGT/Mbp\tTaxon\n')
for each_g in group_list_uniq:
    for_out = '%s\t%s\t%s\t%s\t%s\n' % (each_g, float("{0:.3f}".format(group_to_length_dict[each_g])), group_list_uniq_count[n], group_list_uniq_count_normalized[n], group_2_taxon_dict[each_g])
    HGT_PG_STAT_handle.write(for_out)
    n += 1
HGT_PG_STAT_handle.close()


color_list_according_group_uniq = []
for each in color_list_according_group:
    if color_list_according_group_uniq == []:
        color_list_according_group_uniq.append(each)
    elif each != color_list_according_group_uniq [-1]:
        color_list_according_group_uniq.append(each)


# set xticks fontsize for genome plot
xticks_fontsize_genome = 8
if 25 < len(genome_list_according_group_with_group) <= 50:
    xticks_fontsize_genome = 5
elif 50 < len(genome_list_according_group_with_group) <= 100:
    xticks_fontsize_genome = 4
elif 100 < len(genome_list_according_group_with_group) <= 500:
    xticks_fontsize_genome = 3
elif len(genome_list_according_group_with_group) > 500:
    xticks_fontsize_genome = 2


# set xticks fontsize for group plot
xticks_fontsize_group = 8
if 25 < len(group_list_uniq) <= 50:
    xticks_fontsize_group = 5
elif len(group_list_uniq) > 50:
    xticks_fontsize_group = 3


# set figure size
plt.figure(figsize=(20, 10))
x_range_genome = range(len(num_list_according_group))
x_range_group = range(len(group_list_uniq))

# subplot 1
plt.subplot(221)
plt.bar(x_range_genome, num_list_according_group, alpha=0.5, linewidth=0, color=color_list_according_group)
#plt.xticks([])
plt.xticks(x_range_genome, genome_list_according_group_with_group, rotation=315, fontsize=xticks_fontsize_genome, horizontalalignment='left')
plt.ylabel('Number of HGT')

# subplot 2
plt.subplot(222)
plt.bar(x_range_group, group_list_uniq_count, tick_label=group_id_with_taxon, align='center', alpha=0.5, color=color_list_according_group_uniq)
#plt.xticks([])
plt.xticks(x_range_group, group_id_with_taxon, rotation=315, fontsize=xticks_fontsize_group,horizontalalignment='left')

# subplot 3
plt.subplot(223)
plt.bar(x_range_genome, num_list_according_group_norm, alpha=0.5, linewidth=0, color=color_list_according_group)
#plt.xticks([])
plt.xticks(x_range_genome, genome_list_according_group_with_group, rotation=315, fontsize=xticks_fontsize_genome, horizontalalignment='left')
plt.xlabel('Genome')
plt.ylabel('Number of HGT / Mbp sequences')

# subplot 4
plt.subplot(224)
plt.bar(x_range_group, group_list_uniq_count_normalized, tick_label=group_id_with_taxon, align='center', alpha=0.5, color=color_list_according_group_uniq)
plt.xlabel('Group')
plt.xticks(x_range_group, group_id_with_taxon, rotation=315, fontsize=xticks_fontsize_group,horizontalalignment='left')

# plot layout
#plt.suptitle('Number of HGT detected from each genome/group', fontsize=16)
plt.subplots_adjust(wspace=0.2, top=0.95)

# save plot
plt.tight_layout()
plt.savefig(pwd_candidates_file_ET_validated_STAT_png, dpi=300)


################################################### Get_circlize_plot ##################################################

# for predicted HGTs:
# 1. not at ends
# 2. PG validated

name2id_dict = {}
transfers = []
for each in open(pwd_candidates_file_ET_validated):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        Gene_1 = each_split[0]
        Gene_2 = each_split[1]
        Genome_1_ID = each_split[2]
        Genome_1 = '_'.join(Gene_1.split('_')[:-1])
        Genome_2_ID = each_split[3]
        Genome_2 = '_'.join(Gene_2.split('_')[:-1])
        Identity = each_split[4]
        End_break = each_split[5]
        Direction = each_split[6]
        if Genome_1 not in name2id_dict:
            name2id_dict[Genome_1] = Genome_1_ID
        if Genome_2 not in name2id_dict:
            name2id_dict[Genome_2] = Genome_2_ID
        transfers.append(Direction)

tmp1 = open('tmp1.txt', 'w')
all_group_id = []
for each_t in transfers:
    each_t_split = each_t.split('-->')
    donor = each_t_split[0]
    recipient = each_t_split[1]
    donor_id = name2id_dict[donor].split('_')[0]
    recipient_id = name2id_dict[recipient].split('_')[0]
    if donor_id not in all_group_id:
        all_group_id.append(donor_id)
    if recipient_id not in all_group_id:
        all_group_id.append(recipient_id)
    tmp1.write('%s,%s\n' % (donor_id, recipient_id))
tmp1.close()

os.system('cat tmp1.txt | sort > tmp1_sorted.txt')

current_t = ''
count = 0
tmp2 = open('tmp1_sorted_count.txt', 'w')
for each_t2 in open('tmp1_sorted.txt'):
    each_t2 = each_t2.strip()
    if current_t == '':
        current_t = each_t2
        count += 1
    elif current_t == each_t2:
        count += 1
    elif current_t != each_t2:
        tmp2.write('%s,%s\n' % (current_t, count))
        current_t = each_t2
        count = 1
tmp2.write('%s,%s\n' % (current_t, count))
tmp2.close()

# read in count as dict
transfer_count = {}
for each_3 in open('tmp1_sorted_count.txt'):
    each_3_split = each_3.strip().split(',')
    key = '%s,%s' % (each_3_split[0], each_3_split[1])
    value = each_3_split[2]
    transfer_count[key] = value

all_group_id = sorted(all_group_id)

input_file_name, input_file_ext = os.path.splitext(pwd_candidates_file_ET_validated)
matrix_filename = '%s_matrix.csv' % input_file_name
matrix_file = open(matrix_filename, 'w')
matrix_file.write('\t' + '\t'.join(all_group_id) + '\n')
for each_1 in all_group_id:
    row = [each_1]
    for each_2 in all_group_id:
        current_key = '%s,%s' % (each_2, each_1)
        if current_key not in transfer_count:
            row.append('0')
        else:
            row.append(transfer_count[current_key])
    matrix_file.write('\t'.join(row) + '\n')
matrix_file.close()

# get plot with R
print('Rscript %s -m %s -p %s' % (circos_HGT_R, matrix_filename, pwd_plot_circos))
os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, matrix_filename, pwd_plot_circos))


################################################### remove tmp files ###################################################

# remove tmp files
os.remove('tmp1.txt')
os.remove('tmp1_sorted.txt')
os.remove('tmp1_sorted_count.txt')
os.remove(matrix_filename)

