#!/usr/bin/env python
import os
import shutil
import configparser
import numpy as np
from sys import stdout
from time import sleep
from Bio import SeqIO
from lib.distribution_ploter import plot_identity_list
from lib.get_candidates import get_candidates
from lib.act_ploter import get_gbk_blast_act

############################################## Read in configuration file ##############################################

config = configparser.ConfigParser()
config.read('/Users/songweizhi/Desktop/MetaCHIP_wd/simulated_dataset/original/config.txt')
working_directory = config['FILES_AND_PARAMETERS']['working_directory']

inputs_folder_name = config['FILES_AND_PARAMETERS']['inputs_folder_name']
grouping_file = config['FILES_AND_PARAMETERS']['grouping_file']
cover_cutoff = int(config['FILES_AND_PARAMETERS']['cover_cutoff'])
identity_percentile = int(config['FILES_AND_PARAMETERS']['identity_percentile'])
align_len_cutoff = int(config['FILES_AND_PARAMETERS']['align_len_cutoff'])

blast_results = config['FILES_AND_PARAMETERS']['blast_results']
gbk_file_name = config['FILES_AND_PARAMETERS']['gbk_file_name']

########################################################################################################################


def get_all_identity_list(blast_results, genome_list, alignment_length_cutoff, coverage_cutoff, all_qualified_identities_file):
    # get total match number
    total_match_number = 0
    matches = open(blast_results)
    for match in matches:
        total_match_number += 1

    # get all qualified identities
    matches = open(blast_results)
    out_temp = open(all_qualified_identities_file, 'w')
    all_identities = []
    counted_match = []
    n = 1
    float("{0:.2f}".format(total_match_number/1000))
    for match in matches:
        stdout.write('\r%s x 1000 blast matches detected in total, filtering the %s x 1000th' % (float("{0:.2f}".format(total_match_number/1000)), float("{0:.2f}".format(n/1000))))
        match_split = match.strip().split('\t')
        query = match_split[0]
        subject = match_split[1]
        identity = float(match_split[2])
        align_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        query_split = query.split('_')
        subject_split = subject.split('_')
        query_bin_name = '_'.join(query_split[:-1])
        subject_bin_name = '_'.join(subject_split[:-1])

        coverage_q = float("{0:.2f}".format(float(align_len)*100/float(query_len)))
        coverage_s = float("{0:.2f}".format(float(align_len)*100/float(subject_len)))
        query_name_subject_name = '%s_%s' % (query, subject)
        subject_name_query_name = '%s_%s' % (subject, query)

        # only work on genomes included in grouping file
        if (query_bin_name in genome_list) and (subject_bin_name in genome_list):
            # filter
            if (query_bin_name != subject_bin_name) \
                    and (align_len >= int(alignment_length_cutoff)) \
                    and (coverage_q >= int(coverage_cutoff)) \
                    and (coverage_s >= int(coverage_cutoff)):
                # remove the same match but with swapped query-subject position
                if (query_name_subject_name in counted_match) or (subject_name_query_name in counted_match):
                    pass
                else:
                    all_identities.append(identity)
                    out_temp.write(match)
                    counted_match.append(query_name_subject_name)
            else:
                pass
        else:
            pass
        n += 1

    out_temp.close()
    return all_identities


def do(choice):
    current_group_pair_identities_array = np.array(current_group_pair_identities)
    current_group_pair_identity_cut_off = np.percentile(current_group_pair_identities_array, identity_percentile)
    current_group_pair_identity_cut_off = float("{0:.2f}".format(current_group_pair_identity_cut_off))
    current_group_pair_name_split = current_group_pair_name.split('_')
    current_group_pair_name_swapped = '%s_%s' % (current_group_pair_name_split[1], current_group_pair_name_split[0])
    group_pair_iden_cutoff_dict[current_group_pair_name] = current_group_pair_identity_cut_off
    group_pair_iden_cutoff_dict[current_group_pair_name_swapped] = current_group_pair_identity_cut_off
    if current_group_pair_name == current_group_pair_name_swapped :
        pass
    else:
        group_pair_iden_cutoff_file.write('%s\t%s\n' % (current_group_pair_name, current_group_pair_identity_cut_off))

    if choice in ['s', 'S']:
        pass
    elif choice in ['p', 'P']:
        # check length
        if len(current_group_pair_identities) >= minimum_plot_number:
            if current_group_pair_name == current_group_pair_name_swapped:
                plot_identity_list(current_group_pair_identities, 'None', current_group_pair_name, pwd_iden_distrib_plot_folder)
            else:
                plot_identity_list(current_group_pair_identities, current_group_pair_identity_cut_off, current_group_pair_name, pwd_iden_distrib_plot_folder)
            stdout.write("\rProcessing %dth group-pair: %s, plotting..." % (ploted_group, current_group_pair_name,))
        else:
            unploted_groups.write('%s\t%s\n' % (current_group_pair_name, len(current_group_pair_identities)))
            stdout.write("\rProcessing %dth group-pair: %s, blast match number < %d, skipped" % (ploted_group, current_group_pair_name, minimum_plot_number))


def get_hits_group(input_file_name, output_file_name):
    matches_2 = open(input_file_name)
    output_2_file = open(output_file_name, 'w')
    current_gene = ''
    group_member = []
    for match in matches_2:
        match_split = match.strip().split('\t')
        query_2 = match_split[0]
        target_2 = match_split[1]
        if current_gene == '':
            current_gene = query_2
            group_member.append(target_2)
        else:
            if query_2 == current_gene:
                if target_2 not in group_member:
                    group_member.append(target_2)
                else:
                    pass
            else:
                output_2_file.write('%s\t%s' % (current_gene, '\t'.join(group_member)) + '\n')
                current_gene = query_2
                group_member = []
                group_member.append(target_2)
    output_2_file.write('%s\t%s' % (current_gene, '\t'.join(group_member)) + '\n')
    output_2_file.close()


############################################### Define folder/file name ################################################

print('Define folder/file names and create output folder')

op_prefix = 'output_ip' + str(identity_percentile) + '%_al' + str(align_len_cutoff) + 'bp_c' + str(cover_cutoff) + '%'
op_prefix_iden_0 = 'output_al' + str(align_len_cutoff) + 'bp_c' + str(cover_cutoff) + '%'
op_folder = op_prefix
wd = working_directory

iden_distrib_plot_folder =              'Identity_distribution_images'
qual_idens_file =                       'qualified_identities.txt'
qual_idens_file_gg =                    'qualified_identities_gg.txt'
qual_idens_file_gg_sorted =             'qualified_identities_gg_sorted.txt'
unploted_groups_file =                  'unploted_groups.txt'
qual_idens_with_group_filename =        'qualified_identities_with_group.txt'
qual_idens_with_group_sorted_filename = 'qualified_identities_with_group_sorted.txt'
subjects_in_one_line_filename =         'subjects_in_one_line.txt'
group_pair_iden_cutoff_file_name =      'Identity_cutoff.txt'
op_candidates_with_group_file_name =    'HGT_candidates_with_group.txt'
op_candidates_only_gene_file_name =     'HGT_candidates.txt'
gbk_subset_file =                       'combined_subset.gbk'
op_act_folder_name =                    'ACT_images'

pwd_iden_distrib_plot_folder =      '%s/%s/%s'    % (wd, op_folder, iden_distrib_plot_folder)
pwd_qual_iden_file =                '%s/%s/%s'    % (wd, op_folder, qual_idens_file)
pwd_qual_iden_file_gg =             '%s/%s/%s'    % (wd, op_folder, qual_idens_file_gg)
pwd_qual_iden_file_gg_sorted =      '%s/%s/%s'    % (wd, op_folder, qual_idens_file_gg_sorted)
pwd_unploted_groups_file =          '%s/%s/%s/%s' % (wd, op_folder, iden_distrib_plot_folder, unploted_groups_file)
pwd_qual_idens_with_group =         '%s/%s/%s'    % (wd, op_folder, qual_idens_with_group_filename)
pwd_qual_idens_with_group_sorted =  '%s/%s/%s'    % (wd, op_folder, qual_idens_with_group_sorted_filename)
pwd_subjects_in_one_line =          '%s/%s/%s'    % (wd, op_folder, subjects_in_one_line_filename)
pwd_group_pair_iden_cutoff_file =   '%s/%s/%s'    % (wd, op_folder, group_pair_iden_cutoff_file_name)
pwd_op_candidates_with_group_file = '%s/%s/%s'    % (wd, op_folder, op_candidates_with_group_file_name)
pwd_op_candidates_only_gene_file =  '%s/%s/%s'    % (wd, op_folder, op_candidates_only_gene_file_name)
pwd_gbk_subset_file =               '%s/%s/%s'    % (wd, op_folder, gbk_subset_file)
pwd_op_act_folder =                 '%s/%s/%s'    % (wd, op_folder, op_act_folder_name)
path_to_grouping_file =             '%s/%s/%s'    % (wd, inputs_folder_name, grouping_file)
path_to_blast_results =             '%s/%s/%s'    % (wd, inputs_folder_name, blast_results)
path_to_gbk_file =                  '%s/%s/%s'    % (wd, inputs_folder_name, gbk_file_name)

########################################################################################################################

# forward to working directory
os.chdir(wd)

# create outputs folder
if not os.path.exists(op_folder):
    os.makedirs(op_folder)
else:
    shutil.rmtree(op_folder)
    os.makedirs(op_folder)

# create folder to hold group-group identity distribution plot
os.makedirs(pwd_iden_distrib_plot_folder)


# create genome_group_dict and genome_list
grouping = open(path_to_grouping_file)
genome_name_list = []
name_to_group_number_dict = {}
name_to_group_dict = {}
for each_bin in grouping:
    each_bin_split = each_bin.strip().split(',')
    bin_name = each_bin_split[1]
    bin_group_number = each_bin_split[0]
    bin_group = bin_group_number.split('_')[0]
    genome_name_list.append(bin_name)
    name_to_group_number_dict[bin_name] = bin_group_number
    name_to_group_dict[bin_name] = bin_group

####################################################### Main code ######################################################
sleep(1.5)
print('Plotting overall identity distribution')

# get qualified identities after alignment length and coverage filter (self-match excluded)
all_identities = get_all_identity_list(path_to_blast_results, genome_name_list, align_len_cutoff, cover_cutoff, pwd_qual_iden_file)

# print('debug!!!>_<!!!')

# plot overall identity distribution
all_identities_plot_title = 'All_vs_All'
plot_identity_list(all_identities, 'None', all_identities_plot_title, pwd_iden_distrib_plot_folder)

sleep(1.5)
print('\nPlot overall identity distribution finished')

sleep(1.5)
choice = str(input('Press "P/p" to plot identity distribution for each group,\nPress "S/s" to skip this step:'))

if choice in ['s', 'S']:
    pass
elif choice in ['p', 'P']:
    # Plot identity distribution for each group
    sleep(1.5)
    print('Plotting identity distribution for each group')


# replace query and subject name with query_group-subject_group and sort new file
qualified_identities = open(pwd_qual_iden_file)
qualified_identities_g_g = open(pwd_qual_iden_file_gg, 'w')
for each_identity in qualified_identities:
    each_identity_split = each_identity.strip().split('\t')
    query = each_identity_split[0]
    subject = each_identity_split[1]
    identity = float(each_identity_split[2])
    query_split = query.split('_')
    subject_split = subject.split('_')
    query_genome_name = '_'.join(query_split[:-1])
    subject_genome_name = '_'.join(subject_split[:-1])
    query_group = name_to_group_number_dict[query_genome_name].split('_')[0]
    subject_group = name_to_group_number_dict[subject_genome_name].split('_')[0]
    paired_group_list = [query_group, subject_group]
    # query and subjects name sorted by alphabet order here
    paired_group_list_sorted = sorted(paired_group_list)
    g_g = '%s_%s' % (paired_group_list_sorted[0], paired_group_list_sorted[1])
    qualified_identities_g_g.write('%s\t%s\n' % (g_g, identity))
qualified_identities_g_g.close()

os.system('cat %s | sort > %s' % (pwd_qual_iden_file_gg, pwd_qual_iden_file_gg_sorted))

# get identities for each group pair, plot identity distribution and generate group_pair to identity dict
qualified_identities_g = open(pwd_qual_iden_file_gg_sorted)
unploted_groups = open(pwd_unploted_groups_file, 'w')
current_group_pair_name = ''
current_group_pair_identities = []
group_pair_iden_cutoff_dict = {}
ploted_group = 1
minimum_plot_number = 10
group_pair_iden_cutoff_file = open(pwd_group_pair_iden_cutoff_file, 'w')
for each_identity_g in qualified_identities_g:
    each_identity_g_split = each_identity_g.strip().split('\t')
    group_pair = each_identity_g_split[0]
    identity = float(each_identity_g_split[1])
    group_pair_split = group_pair.split('_')
    group_pair_swapped = '%s_%s' % (group_pair_split[1], group_pair_split[0])
    if current_group_pair_name == '':
        current_group_pair_name = group_pair
        current_group_pair_identities.append(identity)
    else:
        if (group_pair == current_group_pair_name) or (group_pair_swapped == current_group_pair_name):
            current_group_pair_identities.append(identity)
        else:
            # get group_pair to identity dict and identity cut off for defined percentile
            do(choice)
            ploted_group += 1
            if choice in ['s', 'S']:
                pass
            elif choice in ['p', 'P']:
                sleep(0.5)
            # restore current_group_pair_name and current_group_pair_identities for next group pair
            current_group_pair_name = group_pair
            current_group_pair_identities = []
            current_group_pair_identities.append(identity)
# for the last group
do(choice)
group_pair_iden_cutoff_file.close()


sleep(1.5)
if choice in ['p', 'P']:
    print('\nAnalyzing Blast matches to get HGT candidates')
else:
    print('Analyzing Blast matches to get HGT candidates')


# add group information to qualified identities and sort it according to group.
qualified_identities_no_cutoff = open(pwd_qual_iden_file)
qualified_matches_with_group = open(pwd_qual_idens_with_group, 'w')
for qualified_identity in qualified_identities_no_cutoff:
    qualified_identity_split = qualified_identity.strip().split('\t')
    query = qualified_identity_split[0]
    query_split = query.split('_')
    query_bin = '_'.join(query_split[:-1])
    subject = qualified_identity_split[1]
    subject_split = subject.split('_')
    subject_bin = '_'.join(subject_split[:-1])
    identity = float(qualified_identity_split[2])
    file_write = '%s|%s\t%s|%s|%s\n' % (name_to_group_number_dict[query_bin], query, name_to_group_number_dict[subject_bin], subject, str(identity))
    qualified_matches_with_group.write(file_write)
qualified_matches_with_group.close()
# sort qualified_matches_with_group
os.system('cat %s | sort > %s' % (pwd_qual_idens_with_group, pwd_qual_idens_with_group_sorted))

# put all subjects of the same query in one row
get_hits_group(pwd_qual_idens_with_group_sorted, pwd_subjects_in_one_line)

# get HGT candidates
get_candidates(pwd_subjects_in_one_line, pwd_op_candidates_with_group_file, pwd_op_candidates_only_gene_file, group_pair_iden_cutoff_dict)

sleep(1.5)
print('Get HGT candidates finished and exported to %s and %s' % (op_candidates_with_group_file_name, op_candidates_only_gene_file_name))


#################################################### Get ACT images ####################################################

sleep(1.5)
print('Preparing subset of combined gbk file for ACT plotting, be patient!')

# get gene list of all candidates
all_candidates_genes = []
matches = open(pwd_op_candidates_only_gene_file)
for match in matches:
    match_split = match.strip().split('\t')
    for gene in match_split:
        if gene not in all_candidates_genes:
            all_candidates_genes.append(gene)

# print('!!!!!!')
# print(all_candidates_genes)
# print(len(all_candidates_genes))
# print('!!!!!!')

# get subset of combined gbk file
gbk_subset = open(pwd_gbk_subset_file, 'w')
records = SeqIO.parse(path_to_gbk_file, 'genbank')
records_recorded = []





for record in records:
    record_id = record.id
    for gene_f in record.features:
        if 'locus_tag' in gene_f.qualifiers:
            #print(gene_f.qualifiers["locus_tag"])
            for gene_r in all_candidates_genes:
                if gene_r in gene_f.qualifiers["locus_tag"]:
                    if record_id not in records_recorded:
                        SeqIO.write(record, gbk_subset, 'genbank')
                        records_recorded.append(record_id)
                    else:
                        pass

gbk_subset.close()





# create folder to hold ACT output
os.makedirs(pwd_op_act_folder)

sleep(1.5)
print('Get gbk files, run Blast, and plot ACT images')
# plot ACT
gene_gc_content_dict = get_gbk_blast_act(pwd_op_candidates_only_gene_file, pwd_gbk_subset_file, name_to_group_number_dict, pwd_op_act_folder)


sleep(1.5)
print('\nRemove temporary files... ')
# remove temporary files
os.remove(pwd_qual_iden_file)
os.remove(pwd_qual_iden_file_gg)
os.remove(pwd_qual_iden_file_gg_sorted)
os.remove(pwd_qual_idens_with_group)
os.remove(pwd_qual_idens_with_group_sorted)
os.remove(pwd_subjects_in_one_line)
os.remove(pwd_gbk_subset_file)

sleep(1.5)
print('All done, enjoy!')
