import os
import sys
import argparse
import numpy as np
from time import sleep
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime


def get_number_of_group(grouping_file):

    group_list = []
    for each_genome in open(grouping_file):
        each_genome_split = each_genome.strip().split(',')
        group_id = each_genome_split[0]
        if group_id not in group_list:
            group_list.append(group_id)
    number_of_group = len(group_list)

    return number_of_group


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


############################################## input ##############################################

parser = argparse.ArgumentParser()

parser.add_argument('-p', required=True, help='output prefix')
parser.add_argument('-g', required=False, default=None, help='grouping file')
parser.add_argument('-l', required=True, help='grouping level')
parser.add_argument('-cov', required=False, type=int, default=70, help='coverage cutoff')
parser.add_argument('-al', required=False, type=int, default=200, help='alignment length cutoff')
parser.add_argument('-ip', required=False, type=int, default=90, help='identity percentile')
parser.add_argument('-eb', required=False, type=int, default=500, help='the minimal length to be considered as end break')

args = vars(parser.parse_args())
output_prefix = args['p']
grouping_file = args['g']
grouping_level = args['l']
cover_cutoff = args['cov']
align_len_cutoff = args['al']
identity_percentile = args['ip']
ending_match_length = args['eb']


# get path to current script
pwd_get_stat_script = sys.argv[0]
get_stat_script_path, file_name = os.path.split(pwd_get_stat_script)

circos_HGT_R =                      '%s/circos_HGT.R'                     % get_stat_script_path
MetaCHIP_wd =                       '%s_MetaCHIP_wd'                      % output_prefix
MetaCHIP_op_folder =                '%s_HGTs_ip%s_al%sbp_c%s_e%sbp_g%s'   % (output_prefix, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(ending_match_length), get_number_of_group(grouping_file))
candidates_file_name_ET =           'HGT_candidates_PG.txt'
candidates_file_name_ET_validated = 'HGT_candidates_PG_validated.txt'
grouping_id_to_taxon_file_name =    '%s_group_to_taxon_%s.txt'            % (output_prefix, grouping_level)
pwd_candidates_file_ET =            '%s/%s/%s'                            % (MetaCHIP_wd, MetaCHIP_op_folder, candidates_file_name_ET)
pwd_candidates_file_ET_validated =  '%s/%s/%s'                            % (MetaCHIP_wd, MetaCHIP_op_folder, candidates_file_name_ET_validated)
pwd_grouping_id_to_taxon_file =     '%s/%s'                               % (MetaCHIP_wd, grouping_id_to_taxon_file_name)
plot_identity_distribution_BM =     '%s/%s/%s_plot_HGT_identities_BM.png' % (MetaCHIP_wd, MetaCHIP_op_folder, output_prefix)
plot_identity_distribution_PG =     '%s/%s/%s_plot_HGT_identities_PG.png' % (MetaCHIP_wd, MetaCHIP_op_folder, output_prefix)
plot_at_ends_number =               '%s/%s/%s_plot_at_ends_stat.png'      % (MetaCHIP_wd, MetaCHIP_op_folder, output_prefix)
plot_circos =                       '%s/%s/%s_plot_circos_PG.png'         % (MetaCHIP_wd, MetaCHIP_op_folder, output_prefix)

time_format = '[%Y-%m-%d %H:%M:%S] '


############################################## read in prediction results ##############################################

# get identity list
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

print(datetime.now().strftime(time_format) + 'Plot at_ends stat')

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
plt.savefig(plot_at_ends_number, dpi=300)
plt.close()


#################################### plot identity distribution of identified HGTs #####################################

print(datetime.now().strftime(time_format) + 'Plot identity distribution of identified HGTs')

# for predicted HGTs:
# 1. not at ends

# plot Identity distribution of identified HGTs
plot_identity_distribution(identity_list_not_end_BM, plot_identity_distribution_BM)
plot_identity_distribution(identity_list_not_end_PG, plot_identity_distribution_PG)


################################################### Get_circlize_plot ##################################################

print(datetime.now().strftime(time_format) + 'Plot gene flow between groups')

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
#print('Running: Rscript %s -m %s -p %s' % (circos_HGT_R, matrix_filename, circos_plot_name))
os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, matrix_filename, plot_circos))

# remove tmp files
print(datetime.now().strftime(time_format) + 'Deleting temporary files')
os.remove('tmp1.txt')
os.remove('tmp1_sorted.txt')
os.remove('tmp1_sorted_count.txt')
os.remove(matrix_filename)

print(datetime.now().strftime(time_format) + 'Get stat done!')
