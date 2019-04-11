#!/usr/bin/env python

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com or t.thomas@unsw.edu.au

# MetaCHIP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MetaCHIP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import glob
import shutil
import warnings
import argparse
import itertools
import numpy as np
import multiprocessing as mp
from time import sleep
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm
from datetime import datetime
from string import ascii_uppercase
from scipy.stats import gaussian_kde
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from MetaCHIP.MetaCHIP_config import config_dict
# from PIL import Image


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def get_program_path_dict(pwd_cfg_file):
    program_path_dict = {}
    for each in open(pwd_cfg_file):
        each_split = each.strip().split('=')
        program_name = each_split[0]
        program_path = each_split[1]

        # remove space if there are
        if program_name[-1] == ' ':
            program_name = program_name[:-1]
        if program_path[0] == ' ':
            program_path = program_path[1:]

        program_path_dict[program_name] = program_path

    return program_path_dict


def get_group_num_from_grouping_file(grouping_file):
    group_list = set()
    for each_line in open(grouping_file):
        each_line_split = each_line.strip().split(',')
        group_list.add(each_line_split[0])

    group_num = len(group_list)

    return group_num


def get_number_of_group(grouping_file):

    group_list = []
    for each_genome in open(grouping_file):
        each_genome_split = each_genome.strip().split(',')
        group_id = each_genome_split[0]
        if group_id not in group_list:
            group_list.append(group_id)
    number_of_group = len(group_list)

    return number_of_group


def get_group_index_list():
    def iter_all_strings():
        size = 1
        while True:
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)
            size += 1

    group_index_list = []
    for s in iter_all_strings():
        group_index_list.append(s)
        if s == 'ZZ':
            break
    return group_index_list


def index_grouping_file(input_file, output_file):

    output_grouping_with_index = open(output_file, 'w')
    current_group = ''
    current_index = 1
    for each_genome in open(input_file):
        each_genome_split = each_genome.strip().split(',')
        group_id = each_genome_split[0]
        genome_id = each_genome_split[1]
        if current_group == '':
            current_group = group_id
            current_index = 1
            output_grouping_with_index.write('%s_%s,%s\n' % (group_id, current_index, genome_id))
        elif current_group == group_id:
            current_index += 1
            output_grouping_with_index.write('%s_%s,%s\n' % (group_id, current_index, genome_id))
        elif current_group != group_id:
            current_group = group_id
            current_index = 1
            output_grouping_with_index.write('%s_%s,%s\n' % (group_id, current_index, genome_id))


def cluster_2_grouping_file(cluster_file, grouping_file):
    print(os.getcwd())

    t1 = 'cluster_tmp1.txt'
    t1_sorted = 'cluster_tmp1_sorted.txt'

    t1_handle = open(t1, 'w')
    for each in open(cluster_file):
        if not each.startswith(','):
            genome_id = each.strip().split(',')[0]
            cluster_id = each.strip().split(',')[1]
            t1_handle.write('%s,%s\n' % (cluster_id, genome_id))
    t1_handle.close()

    os.system('cat %s | sort > %s' % (t1, t1_sorted))
    group_index_list = get_group_index_list()
    grouping_file_handle = open(grouping_file, 'w')

    current_cluster_name = ''
    group_index_no = 0
    n = 1
    for each in open(t1_sorted):
        cluster_name = each.strip().split(',')[0]
        genome_name = each.strip().split(',')[1]

        if current_cluster_name == '':
            current_cluster_name = cluster_name
            grouping_file_handle.write('%s_%s,%s\n' % (group_index_list[group_index_no], n, genome_name))
            n += 1
        elif current_cluster_name == cluster_name:
            grouping_file_handle.write('%s_%s,%s\n' % (group_index_list[group_index_no], n, genome_name))
            n += 1
        elif current_cluster_name != cluster_name:
            current_cluster_name = cluster_name
            group_index_no += 1
            n = 1
            grouping_file_handle.write('%s_%s,%s\n' % (group_index_list[group_index_no], n, genome_name))
            n += 1

    os.remove(t1)
    os.remove(t1_sorted)


def uniq_list(input_list):
    output_list = []
    for each_element in input_list:
        if each_element not in output_list:
            output_list.append(each_element)
    return output_list


def get_qualigied_blast_hits(pwd_blast_results, align_len_cutoff, cover_cutoff, genome_name_list, pwd_qual_iden_file):

    out_temp = open(pwd_qual_iden_file, 'w')
    for match in open(pwd_blast_results):
        match_split = match.strip().split('\t')
        query = match_split[0]
        subject = match_split[1]
        align_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        query_bin_name = '_'.join(query.split('_')[:-1])
        subject_bin_name = '_'.join(subject.split('_')[:-1])
        coverage_q = float(align_len) * 100 / float(query_len)
        coverage_s = float(align_len) * 100 / float(subject_len)

        # first filter with alignment length
        if align_len >= int(align_len_cutoff):

            # then remove within genome hits
            if query_bin_name != subject_bin_name:

                # then coverage cutoff
                if (coverage_q >= int(cover_cutoff)) and (coverage_s >= int(cover_cutoff)):

                    # only work on genomes with clear taxonomic classification
                    if (query_bin_name in genome_name_list) and (subject_bin_name in genome_name_list):
                        out_temp.write(match)
    out_temp.close()


def plot_identity_list(identity_list, identity_cut_off, title, output_foler):
    identity_list = sorted(identity_list)

    # get statistics
    match_number = len(identity_list)
    average_iden = float(np.average(identity_list))
    average_iden = float("{0:.2f}".format(average_iden))
    max_match = float(np.max(identity_list))
    min_match = float(np.min(identity_list))

    # get hist plot
    num_bins = 50
    plt.hist(identity_list, num_bins, alpha=0.1, normed=1, facecolor='blue')  # normed = 1 normalized to 1, that is probablity
    plt.title('Group: %s' % title)
    plt.xlabel('Identity')
    plt.ylabel('Probability')
    plt.subplots_adjust(left=0.15)

    # get fit line
    density = gaussian_kde(identity_list)
    x_axis = np.linspace(min_match - 5, max_match + 5, 200)
    density.covariance_factor = lambda: 0.3
    density._compute_covariance()
    plt.plot(x_axis, density(x_axis))

    # add text
    x_min = plt.xlim()[0]  # get the x-axes minimum value
    x_max = plt.xlim()[1]  # get the x-axes maximum value
    y_min = plt.ylim()[0]  # get the y-axes minimum value
    y_max = plt.ylim()[1]  # get the y-axes maximum value

    # set text position
    text_x = x_min + (x_max - x_min)/5 * 3.8
    text_y_total = y_min + (y_max - y_min) / 5 * 4.4
    text_y_min = y_min + (y_max - y_min) / 5 * 4.1
    text_y_max = y_min + (y_max - y_min) / 5 * 3.8
    text_y_average = y_min + (y_max - y_min) / 5 * 3.5
    text_y_cutoff = y_min + (y_max - y_min) / 5 * 3.2

    # plot text
    plt.text(text_x, text_y_total, 'Total: %s' % match_number)
    plt.text(text_x, text_y_min, 'Min: %s' % min_match)
    plt.text(text_x, text_y_max, 'Max: %s' % max_match)
    plt.text(text_x, text_y_average, 'Mean: %s' % average_iden)
    plt.text(text_x, text_y_cutoff, 'Cutoff: %s' % identity_cut_off)
    if identity_cut_off != 'None':
        plt.annotate(' ',
                     xy=(identity_cut_off, 0),
                     xytext=(identity_cut_off, density(identity_cut_off)),
                     arrowprops=dict(width=0.5,
                                headwidth=0.5,
                                facecolor='red',
                                edgecolor='red',
                                shrink=0.02))
    # Get plot
    plt.savefig('%s/%s.png' % (output_foler, title), dpi = 300)
    plt.close()


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


def get_candidates(targets_group_file, gene_with_g_file_name, gene_only_name_file_name, group_pair_iden_cutoff_dict):

    output_1 = open(gene_with_g_file_name, 'w')
    output_2 = open(gene_only_name_file_name, 'w')

    for group in open(targets_group_file):
        group_split = group.strip().split('\t')
        query = group_split[0]
        query_split = query.split('|')
        query_gene_name = query_split[1]
        query_sg = query_split[0]
        query_g = query_sg.split('_')[0]
        subjects_list = group_split[1:]

        # if only one non-self subject was found, no matter which group it comes from, ignored
        if len(subjects_list) == 1:
            pass

        # if more than 1 non-self subject was found:
        elif len(subjects_list) > 1:
            # get the number of subjects from self-group and non_self_group
            self_group_subject_list = []
            non_self_group_subject_list = []
            #print(non_self_group_subject_list)
            for each_subject in subjects_list:

                #print(each_subject)
                #each_subject_g = each_subject[0]
                each_subject_g = each_subject.split('|')[0].split('_')[0]
                if each_subject_g == query_g:
                    self_group_subject_list.append(each_subject)
                else:
                    non_self_group_subject_list.append(each_subject)


            # if only the self-match was found in self-group, all matched from other groups, if any, will be ignored
            if len(self_group_subject_list) == 0:
                pass

            # if no non-self-group subjects was found, ignored
            elif (len(self_group_subject_list) > 0) and (len(non_self_group_subject_list) == 0):
                pass

            # if both non-self self-group subjects and non-self-group subject exist:
            elif (len(self_group_subject_list) > 0) and (len(non_self_group_subject_list) > 0):
                # get the number the groups
                non_self_group_subject_list_uniq = []
                for each_g in non_self_group_subject_list:
                    each_g_group = each_g.split('|')[0].split('_')[0]
                    if each_g_group not in non_self_group_subject_list_uniq:
                        non_self_group_subject_list_uniq.append(each_g_group)

                # if all non-self-group subjects come from the same group
                if len(non_self_group_subject_list_uniq) == 1:

                    # get the maximum and average identity from self-group
                    sg_maximum = 0
                    sg_sum = 0
                    sg_subject_number = 0
                    for each_sg_subject in self_group_subject_list:
                        each_sg_subject_iden = float(each_sg_subject.split('|')[2])
                        if each_sg_subject_iden > sg_maximum:
                            sg_maximum = each_sg_subject_iden
                        sg_sum += each_sg_subject_iden
                        sg_subject_number += 1
                    sg_average = sg_sum/sg_subject_number

                    # get the maximum and average identity from non-self-group
                    nsg_maximum = 0
                    nsg_maximum_gene = ''
                    nsg_sum = 0
                    nsg_subject_number = 0
                    for each_nsg_subject in non_self_group_subject_list:
                        each_nsg_subject_iden = float(each_nsg_subject.split('|')[2])
                        if each_nsg_subject_iden > nsg_maximum:
                            nsg_maximum = each_nsg_subject_iden
                            nsg_maximum_gene = each_nsg_subject
                        nsg_sum += each_nsg_subject_iden
                        nsg_subject_number += 1
                    nsg_average = nsg_sum/nsg_subject_number

                    # if the average non-self-group identity > average self-group identity,
                    # Subject with maximum identity from this group will be considered as a HGT donor.
                    if nsg_average > sg_average:
                        # filter with obtained identity cut-off:
                        candidate_g = nsg_maximum_gene.split('|')[0].split('_')[0]
                        candidate_iden = float(nsg_maximum_gene.split('|')[2])
                        qg_sg = '%s_%s' % (query_g, candidate_g)
                        qg_sg_iden_cutoff = group_pair_iden_cutoff_dict[qg_sg]
                        if candidate_iden >= qg_sg_iden_cutoff:
                            output_1.write('%s\t%s\n' % (query, nsg_maximum_gene))
                            output_2.write('%s\t%s\n' % (query_gene_name, nsg_maximum_gene.split('|')[1]))

                # if non-self-group subjects come from different groups
                elif len(non_self_group_subject_list_uniq) > 1:
                    # get average/maximum for self-group
                    sg_maximum = 0
                    sg_sum = 0
                    sg_subject_number = 0
                    for each_sg_subject in self_group_subject_list:
                        each_sg_subject_iden = float(each_sg_subject.split('|')[2])
                        if each_sg_subject_iden > sg_maximum:
                            sg_maximum = each_sg_subject_iden
                        sg_sum += each_sg_subject_iden
                        sg_subject_number += 1
                    sg_average = sg_sum/sg_subject_number

                    # get average/maximum for each non-self-group
                    nsg_average_dict = {}
                    nsg_maximum_dict = {}
                    nsg_maximum_gene_name_dict = {}
                    for each_nsg in non_self_group_subject_list_uniq:
                        nsg_maximum = 0
                        nsg_maximum_gene = ''
                        nsg_sum = 0
                        nsg_subject_number = 0
                        for each_nsg_subject in non_self_group_subject_list:
                            #print(non_self_group_subject_list)
                            each_nsg_subject_iden = float(each_nsg_subject.split('|')[2])
                            #print(each_nsg_subject)
                            #print(each_nsg_subject[0])
                            each_nsg_subject_group = each_nsg_subject.split('|')[0].split('_')[0]
                            if each_nsg_subject_group == each_nsg:
                                if each_nsg_subject_iden > nsg_maximum:
                                    nsg_maximum = each_nsg_subject_iden
                                    nsg_maximum_gene = each_nsg_subject
                                nsg_sum += each_nsg_subject_iden
                                nsg_subject_number += 1
                        nsg_average = nsg_sum / nsg_subject_number
                        nsg_average_dict[each_nsg] = nsg_average
                        nsg_maximum_dict[each_nsg] = nsg_maximum
                        nsg_maximum_gene_name_dict[each_nsg] = nsg_maximum_gene

                    # get the group with maximum average group identity
                    maximum_average = sg_average
                    maximum_average_g = group[0]
                    for each_g in nsg_average_dict:
                        if nsg_average_dict[each_g] > maximum_average:
                            maximum_average = nsg_average_dict[each_g]
                            maximum_average_g = each_g

                    # if the maximum average identity group is the self-group, ignored
                    if maximum_average_g == group[0]:
                        pass

                    # if self-group average identity is not the maximum,
                    # Group with maximum average identity will be considered as the candidate donor group,
                    # Subject with maximum identity from the candidate donor group will be considered as a HGT donor.
                    elif maximum_average_g != group[0]:
                        # filter with obtained identity cut-off:
                        candidate_g = nsg_maximum_gene_name_dict[maximum_average_g].split('|')[0].split('_')[0]
                        candidate_iden = float(nsg_maximum_gene_name_dict[maximum_average_g].split('|')[2])
                        qg_sg = '%s_%s' % (query_g, candidate_g)
                        qg_sg_iden_cutoff = group_pair_iden_cutoff_dict[qg_sg]
                        if candidate_iden >= qg_sg_iden_cutoff:
                            output_1.write('%s\t%s\n' % (query, nsg_maximum_gene_name_dict[maximum_average_g]))
                            output_2.write('%s\t%s\n' % (query_gene_name, nsg_maximum_gene_name_dict[maximum_average_g].split('|')[1]))
    output_1.close()
    output_2.close()


def check_match_direction(blast_hit_splitted):
    query_start = int(blast_hit_splitted[6])
    query_end = int(blast_hit_splitted[7])
    subject_start = int(blast_hit_splitted[8])
    subject_end = int(blast_hit_splitted[9])
    query_direction = query_end - query_start
    subject_direction = subject_end - subject_start

    same_match_direction = True
    if ((query_direction > 0) and (subject_direction < 0)) or ((query_direction < 0) and (subject_direction > 0)):
        same_match_direction = False

    return same_match_direction


def check_full_lenght_and_end_match(qualified_ctg_match_list, identity_cutoff):

    ######################################## check full length match ########################################

    query_len = int(qualified_ctg_match_list[0][12])
    subject_len = int(qualified_ctg_match_list[0][13])

    # get position list of matched regions
    query_matched_region_list = []
    subject_matched_region_list = []
    for ctg_match in qualified_ctg_match_list:
        query_start = int(ctg_match[6])
        query_end = int(ctg_match[7])
        subject_start = int(ctg_match[8])
        subject_end = int(ctg_match[9])
        query_matched_region = sorted([query_start, query_end])
        subject_matched_region = sorted([subject_start, subject_end])
        query_matched_region_list.append(query_matched_region)
        subject_matched_region_list.append(subject_matched_region)

    # get total length of query matched regions
    query_matched_len_total = 0
    current_query_end = 0
    for query_matched in sorted(query_matched_region_list):

        if query_matched_len_total == 0:
            query_matched_len_total = query_matched[1] - query_matched[0] + 1
            current_query_end = query_matched[1]

        elif query_matched[0] > current_query_end:
            query_matched_len_total += query_matched[1] - query_matched[0] + 1

        elif query_matched[0] < current_query_end:
            if query_matched[1] > current_query_end:
                query_matched_len_total += query_matched[1] - current_query_end

        elif query_matched[0] == current_query_end:
            query_matched_len_total += query_matched[1] - current_query_end

    # get total length of subject matched regions
    subject_matched_len_total = 0
    current_subject_end = 0
    for subject_matched in sorted(subject_matched_region_list):

        if subject_matched_len_total == 0:
            subject_matched_len_total = subject_matched[1] - subject_matched[0] + 1
            current_subject_end = subject_matched[1]

        elif subject_matched[0] > current_subject_end:
            subject_matched_len_total += subject_matched[1] - subject_matched[0] + 1

        elif subject_matched[0] < current_subject_end:
            if subject_matched[1] > current_subject_end:
                subject_matched_len_total += subject_matched[1] - current_subject_end

        elif subject_matched[0] == current_subject_end:
            subject_matched_len_total += subject_matched[1] - current_subject_end


    # get total coverage for query and subject
    query_cov_total = query_matched_len_total/query_len
    subject_cov_total = subject_matched_len_total/subject_len

    # get match category
    match_category = 'normal'
    best_hit_end_gap_len = 200
    gap_cutoff_for_concatenating = 300
    # full length match: coverage cutoff 95%
    if (query_cov_total >= 0.95) or (subject_cov_total >= 0.95):
        match_category = 'full_length_match'


    ######################################## check end match ########################################

    else:

        # read in best hit information
        best_hit = qualified_ctg_match_list[0]
        best_hit_identity = float(best_hit[2])
        best_hit_query_start = int(best_hit[6])
        best_hit_query_end = int(best_hit[7])
        best_hit_subject_start = int(best_hit[8])
        best_hit_subject_end = int(best_hit[9])
        query_len = int(best_hit[12])
        subject_len = int(best_hit[13])
        best_hit_same_direction = check_match_direction(best_hit)

        # concatenate continuously matched blocks with gap less than 200bp
        matched_block_query_start = best_hit_query_start
        matched_block_query_end = best_hit_query_end
        matched_block_subject_start = best_hit_subject_start
        matched_block_subject_end = best_hit_subject_end

        if best_hit_identity >= identity_cutoff:

            for matched_block in qualified_ctg_match_list[1:]:

                current_block_identity = float(matched_block[2])
                current_block_direction = check_match_direction(matched_block)

                # if identity difference <= 1 and has same match direction with the best hit
                if (-6 <= (best_hit_identity - current_block_identity) <= 6) and (current_block_direction == best_hit_same_direction):

                    current_query_start = int(matched_block[6])
                    current_query_end = int(matched_block[7])
                    current_subject_start = int(matched_block[8])
                    current_subject_end = int(matched_block[9])

                    if best_hit_same_direction is True:

                        # situation 1
                        if ((current_query_start >= matched_block_query_start) and (current_query_end <= matched_block_query_end)) and ((current_subject_start >= matched_block_subject_start) and (current_subject_end <= matched_block_subject_end)):
                            pass  # do nothing

                        # situation 2
                        if ((current_query_start > matched_block_query_start) and (current_query_end > matched_block_query_end)) and ((current_subject_start > matched_block_subject_start) and (current_subject_end > matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (current_query_start - matched_block_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (current_subject_start - matched_block_subject_end) <= gap_cutoff_for_concatenating):
                            matched_block_query_end = current_query_end
                            matched_block_subject_end = current_subject_end

                        # situation 3
                        if ((current_query_start < matched_block_query_start) and (current_query_end < matched_block_query_end)) and ((current_subject_start < matched_block_subject_start) and (current_subject_end < matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (matched_block_query_start - current_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (matched_block_subject_start - current_subject_end) <= gap_cutoff_for_concatenating):
                            matched_block_query_start = current_query_start
                            matched_block_subject_start = current_subject_start

                    if best_hit_same_direction is False:

                        # situation 1
                        if ((current_query_start >= matched_block_query_start) and (current_query_end <= matched_block_query_end)) and ((current_subject_start <= matched_block_subject_start) and (current_subject_end >= matched_block_subject_end)):
                            pass  # do nothing

                        # situation 2
                        if ((current_query_start > matched_block_query_start) and (current_query_end > matched_block_query_end)) and ((current_subject_start < matched_block_subject_start) and (current_subject_end < matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (current_query_start - matched_block_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (matched_block_subject_end - current_subject_start) <= gap_cutoff_for_concatenating):
                            matched_block_query_end = current_query_end
                            matched_block_subject_end = current_subject_end

                        # situation 3
                        if ((current_query_start < matched_block_query_start) and (current_query_end < matched_block_query_end)) and ((current_subject_start > matched_block_subject_start) and (current_subject_end > matched_block_subject_end)) and (-gap_cutoff_for_concatenating <= (matched_block_query_start - current_query_end) <= gap_cutoff_for_concatenating) and (-gap_cutoff_for_concatenating <= (current_subject_end - matched_block_subject_start) <= gap_cutoff_for_concatenating):
                            matched_block_query_start = current_query_start
                            matched_block_subject_start = current_subject_start


            ######################################## check end_match ########################################

            # situation 1
            if (best_hit_same_direction is True) and (query_len - matched_block_query_end <= best_hit_end_gap_len) and (matched_block_subject_start <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 2
            elif (best_hit_same_direction is True) and (matched_block_query_start <= best_hit_end_gap_len) and (subject_len - matched_block_subject_end <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 3
            elif (best_hit_same_direction is False) and (query_len - matched_block_query_end <= best_hit_end_gap_len) and (subject_len - matched_block_subject_start <= best_hit_end_gap_len):
                match_category = 'end_match'

            # situation 4
            elif (best_hit_same_direction is False) and (matched_block_query_start <= best_hit_end_gap_len) and (matched_block_subject_end <= best_hit_end_gap_len):
                match_category = 'end_match'

    return match_category


def set_contig_track_features(gene_contig, name_group_dict, candidate_list, HGT_iden, feature_set):
    # add features to feature set
    for feature in gene_contig.features:
        if feature.type == "CDS":

            # define label color
            if feature.qualifiers['locus_tag'][0] in candidate_list:
                label_color = colors.blue
                label_size = 16
                # add identity to gene is
                feature.qualifiers['locus_tag'][0] = '%s (%s)' % (feature.qualifiers['locus_tag'][0], HGT_iden)
            else:
                label_color = colors.black
                label_size = 10

            # change gene name
            bin_name_gbk_split = feature.qualifiers['locus_tag'][0].split('_')
            bin_name_gbk = '_'.join(bin_name_gbk_split[:-1])
            feature.qualifiers['locus_tag'][0] = '%s_%s' % (name_group_dict[bin_name_gbk], bin_name_gbk_split[-1])

            # strands
            color = None
            label_angle = 0
            if feature.location.strand == 1:
                label_angle = 45
                color = colors.lightblue
            elif feature.location.strand == -1:
                label_angle = -225
                color = colors.lightgreen
            # add feature
            feature_set.add_feature(feature,
                                    color=color,
                                    label=True,
                                    sigil='ARROW',
                                    arrowshaft_height=0.5,
                                    arrowhead_length=0.4,
                                    label_color=label_color,
                                    label_size=label_size,
                                    label_angle=label_angle,
                                    label_position="middle")


def get_flanking_region(input_gbk_file, HGT_candidate, flanking_length):

    wd, gbk_file = os.path.split(input_gbk_file)
    new_gbk_file = '%s/%s_%sbp_temp.gbk' % (wd, HGT_candidate, flanking_length)
    new_gbk_final_file = '%s/%s_%sbp.gbk' % (wd, HGT_candidate, flanking_length)
    new_fasta_final_file = '%s/%s_%sbp.fasta' % (wd, HGT_candidate, flanking_length)

    # get flanking range of candidate
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_start = 0
    new_end = 0
    contig_length = 0
    for record in input_gbk:
        contig_length = len(record.seq)
        for gene in record.features:
            # get contig length
            if gene.type == 'source':
                pass
            # get new start and end points
            elif 'locus_tag' in gene.qualifiers:
                if gene.qualifiers['locus_tag'][0] == HGT_candidate:
                    # get new start point
                    new_start = gene.location.start - flanking_length
                    if new_start < 0:
                        new_start = 0
                    # get new end point
                    new_end = gene.location.end + flanking_length
                    if new_end > contig_length:
                        new_end = contig_length

    # get genes within flanking region
    keep_gene_list = []
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    for record in input_gbk:
        for gene in record.features:
            if 'locus_tag' in gene.qualifiers:
                if (gene.location.start < new_start) and (gene.location.end >= new_start):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                    new_start = gene.location.start
                elif (gene.location.start > new_start) and (gene.location.end < new_end):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                elif (gene.location.start <= new_end) and (gene.location.end > new_end):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                    new_end = gene.location.end

    # remove genes not in flanking region from gbk file
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_gbk = open(new_gbk_file, 'w')
    for record in input_gbk:
        new_record_features = []
        for gene in record.features:
            if gene.type == 'source':
                new_record_features.append(gene)
            elif 'locus_tag' in gene.qualifiers:
                if gene.qualifiers['locus_tag'][0] in keep_gene_list:
                    new_record_features.append(gene)
        record.features = new_record_features
        SeqIO.write(record, new_gbk, 'genbank')
    new_gbk.close()

    # remove sequences not in flanking region
    new_gbk_full_length = SeqIO.parse(new_gbk_file, "genbank")
    new_gbk_final = open(new_gbk_final_file, 'w')
    new_fasta_final = open(new_fasta_final_file, 'w')
    for record in new_gbk_full_length:
        # get new sequence
        new_seq = record.seq[new_start:new_end]
        new_contig_length = len(new_seq)
        new_record = SeqRecord(new_seq,
                               id=record.id,
                               name=record.name,
                               description=record.description,
                               annotations=record.annotations)

        # get new location
        new_record_features_2 = []
        for gene in record.features:
            if gene.type == 'source':
                gene_location_new = ''
                if gene.location.strand == 1:
                    gene_location_new = FeatureLocation(0, new_contig_length, strand=+1)
                if gene.location.strand == -1:
                    gene_location_new = FeatureLocation(0, new_contig_length, strand=-1)
                gene.location = gene_location_new
                new_record_features_2.append(gene)
            elif 'locus_tag' in gene.qualifiers:
                gene_location_new = ''
                if gene.location.strand == 1:
                    if gene.location.start - new_start < 0:
                        gene_location_new = FeatureLocation(0,
                                                            gene.location.end - new_start,
                                                            strand=+1)
                    else:
                        gene_location_new = FeatureLocation(gene.location.start - new_start,
                                                            gene.location.end - new_start,
                                                            strand=+1)
                if gene.location.strand == -1:
                    if gene.location.start - new_start < 0:
                        gene_location_new = FeatureLocation(0,
                                                            gene.location.end - new_start,
                                                            strand=-1)
                    else:
                        gene_location_new = FeatureLocation(gene.location.start - new_start,
                                                            gene.location.end - new_start,
                                                            strand=-1)
                gene.location = gene_location_new
                new_record_features_2.append(gene)
        new_record.features = new_record_features_2
        SeqIO.write(new_record, new_gbk_final, 'genbank')
        SeqIO.write(new_record, new_fasta_final, 'fasta')

    new_gbk_final.close()
    new_fasta_final.close()
    os.system('rm %s' % new_gbk_file)


def export_dna_record(gene_seq, gene_id, gene_description, pwd_output_file):
    output_handle = open(pwd_output_file, 'w')
    seq_object = Seq(str(gene_seq), IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')
    output_handle.close()


def get_gbk_blast_act2(arguments_list):

    match = arguments_list[0]
    pwd_gbk_folder = arguments_list[1]
    flanking_length = arguments_list[2]
    aln_len_cutoff = arguments_list[3]
    name_to_group_number_dict = arguments_list[4]
    path_to_output_act_folder = arguments_list[5]
    pwd_normal_plot_folder = arguments_list[6]
    pwd_at_ends_plot_folder = arguments_list[7]
    pwd_full_contig_match_plot_folder = arguments_list[8]
    pwd_blastn_exe = arguments_list[9]
    keep_temp = arguments_list[10]
    candidates_2_contig_match_category_dict = arguments_list[11]
    end_match_iden_cutoff = arguments_list[12]
    No_Eb_Check = arguments_list[13]
    flk_plot_fmt = 'SVG'

    genes = match.strip().split('\t')[:-1]
    current_HGT_iden = float("{0:.1f}".format(float(match.strip().split('\t')[-1])))
    folder_name = '___'.join(genes)
    os.mkdir('%s/%s' % (path_to_output_act_folder, folder_name))

    gene_1 = genes[0]
    gene_2 = genes[1]
    genome_1 = '_'.join(gene_1.split('_')[:-1])
    genome_2 = '_'.join(gene_2.split('_')[:-1])
    pwd_genome_1_gbk = '%s/%s.gbk' % (pwd_gbk_folder, genome_1)
    pwd_genome_2_gbk = '%s/%s.gbk' % (pwd_gbk_folder, genome_2)

    dict_value_list = []
    # Extract gbk and fasta files for gene 1
    for genome_1_record in SeqIO.parse(pwd_genome_1_gbk, 'genbank'):
        for gene_1_f in genome_1_record.features:
            if 'locus_tag' in gene_1_f.qualifiers:
                if gene_1 in gene_1_f.qualifiers["locus_tag"]:
                    dict_value_list.append([gene_1, int(gene_1_f.location.start), int(gene_1_f.location.end), gene_1_f.location.strand, len(genome_1_record.seq)])
                    pwd_gene_1_gbk_file = '%s/%s/%s.gbk' % (path_to_output_act_folder, folder_name, gene_1)
                    pwd_gene_1_fasta_file = '%s/%s/%s.fasta' % (path_to_output_act_folder, folder_name, gene_1)
                    SeqIO.write(genome_1_record, pwd_gene_1_gbk_file, 'genbank')
                    SeqIO.write(genome_1_record, pwd_gene_1_fasta_file, 'fasta')
                    # get flanking regions
                    get_flanking_region(pwd_gene_1_gbk_file, gene_1, flanking_length)

    # Extract gbk and fasta files for gene 2
    for genome_2_record in SeqIO.parse(pwd_genome_2_gbk, 'genbank'):
        for gene_2_f in genome_2_record.features:
            if 'locus_tag' in gene_2_f.qualifiers:
                if gene_2 in gene_2_f.qualifiers["locus_tag"]:
                    dict_value_list.append([gene_2, int(gene_2_f.location.start), int(gene_2_f.location.end), gene_2_f.location.strand, len(genome_2_record.seq)])
                    pwd_gene_2_gbk_file = '%s/%s/%s.gbk' % (path_to_output_act_folder, folder_name, gene_2)
                    pwd_gene_2_fasta_file = '%s/%s/%s.fasta' % (path_to_output_act_folder, folder_name, gene_2)
                    SeqIO.write(genome_2_record, pwd_gene_2_gbk_file, 'genbank')
                    SeqIO.write(genome_2_record, pwd_gene_2_fasta_file, 'fasta')
                    # get flanking regions
                    get_flanking_region(pwd_gene_2_gbk_file, gene_2, flanking_length)

    # Run Blast
    prefix_c =              '%s/%s'                 % (path_to_output_act_folder, folder_name)
    query_c =               '%s/%s_%sbp.fasta'      % (prefix_c, gene_1, flanking_length)
    query_c_full_len =      '%s/%s.fasta'           % (prefix_c, gene_1)
    subject_c =             '%s/%s_%sbp.fasta'      % (prefix_c, gene_2, flanking_length)
    subject_c_full_len =    '%s/%s.fasta'           % (prefix_c, gene_2)
    output_c =              '%s/%s.txt'             % (prefix_c, folder_name)
    output_c_full_len =     '%s/%s_full_length.txt' % (prefix_c, folder_name)

    parameters_c_n =          '-evalue 1e-5 -outfmt 6 -task blastn'
    parameters_c_n_full_len = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
    command_blast =           '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_c, subject_c, output_c, parameters_c_n)
    command_blast_full_len =  '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_c_full_len, subject_c_full_len, output_c_full_len, parameters_c_n_full_len)

    os.system(command_blast)

    if No_Eb_Check is False:
        os.system(command_blast_full_len)

    ############################## check whether full length or end match ##############################

    # get match category
    if No_Eb_Check is True:
        match_category = 'normal'
    else:
        # get qualified_ctg_match_list
        min_ctg_match_aln_len = 100
        qualified_ctg_match_list = []
        for blast_hit in open(output_c_full_len):
            blast_hit_split = blast_hit.strip().split('\t')
            align_len = int(blast_hit_split[3])
            if align_len >= min_ctg_match_aln_len:
                qualified_ctg_match_list.append(blast_hit_split)

        if len(qualified_ctg_match_list) == 0:
            match_category = 'normal'
        else:
            match_category = check_full_lenght_and_end_match(qualified_ctg_match_list, end_match_iden_cutoff)

    candidates_2_contig_match_category_dict[folder_name] = match_category


    ############################## prepare for flanking plot ##############################

    # read in gbk files
    matche_pair_list = []
    for each_gene in genes:
        path_to_gbk_file = '%s/%s/%s_%sbp.gbk' % (path_to_output_act_folder, folder_name, each_gene, flanking_length)
        gene_contig = SeqIO.read(path_to_gbk_file, "genbank")
        matche_pair_list.append(gene_contig)
    bin_record_list = []
    bin_record_list.append(matche_pair_list)

    # get the distance of the gene to contig ends
    gene_1_left_len = dict_value_list[0][1]
    gene_1_right_len = dict_value_list[0][4] - dict_value_list[0][2]
    gene_2_left_len = dict_value_list[1][1]
    gene_2_right_len = dict_value_list[1][4] - dict_value_list[1][2]

    # create an empty diagram
    diagram = GenomeDiagram.Diagram()

    # add tracks to diagram
    max_len = 0
    for gene1_contig, gene2_contig in bin_record_list:
        # set diagram track length
        max_len = max(max_len, len(gene1_contig), len(gene2_contig))

        # add gene content track for gene1_contig
        contig_1_gene_content_track = diagram.new_track(1,
                                                        name='%s (left %sbp, right %sbp)' % (gene1_contig.name, gene_1_left_len, gene_1_right_len),
                                                        greytrack=True,
                                                        greytrack_labels=1,
                                                        greytrack_font='Helvetica',
                                                        greytrack_fontsize=12,
                                                        height=0.35,
                                                        start=0,
                                                        end=len(gene1_contig),
                                                        scale=True,
                                                        scale_fontsize=6,
                                                        scale_ticks=1,
                                                        scale_smalltick_interval=10000,
                                                        scale_largetick_interval=10000)

        # add gene content track for gene2_contig
        contig_2_gene_content_track = diagram.new_track(1,
                                                        name='%s (left %sbp, right %sbp)' % (gene2_contig.name, gene_2_left_len, gene_2_right_len),
                                                        greytrack=True,
                                                        greytrack_labels=1,
                                                        greytrack_font='Helvetica',
                                                        greytrack_fontsize=12,
                                                        height=0.35,
                                                        start=0,
                                                        end=len(gene2_contig),
                                                        scale=True,
                                                        scale_fontsize=6,
                                                        scale_ticks=1,
                                                        scale_smalltick_interval=10000,
                                                        scale_largetick_interval=10000)

        # add blank feature/graph sets to each track
        feature_sets_1 = contig_1_gene_content_track.new_set(type='feature')
        feature_sets_2 = contig_2_gene_content_track.new_set(type='feature')

        # add gene features to 2 blank feature sets
        set_contig_track_features(gene1_contig, name_to_group_number_dict, genes, current_HGT_iden, feature_sets_1)
        set_contig_track_features(gene2_contig, name_to_group_number_dict, genes, current_HGT_iden, feature_sets_2)

        ####################################### add crosslink from blast results #######################################

        path_to_blast_result = '%s/%s/%s.txt' % (path_to_output_act_folder, folder_name, folder_name)
        blast_results = open(path_to_blast_result)

        # parse blast results
        for each_line in blast_results:
            each_line_split = each_line.split('\t')
            query = each_line_split[0]
            identity = float(each_line_split[2])
            alignment_len = int(each_line_split[3])
            query_start = int(each_line_split[6])
            query_end = int(each_line_split[7])
            target_start = int(each_line_split[8])
            target_end = int(each_line_split[9])

            # use color to reflect identity
            color = colors.linearlyInterpolatedColor(colors.white, colors.red, 50, 100, identity)

            # only focus on matches longer than 100 bp
            if alignment_len >= 200:

                # if query is contig_1
                if query == gene1_contig.name:
                    link = CrossLink((contig_1_gene_content_track, query_start, query_end),
                                     (contig_2_gene_content_track, target_start, target_end),
                                     color=color,
                                     border=color,
                                     flip=False)
                    diagram.cross_track_links.append(link)

                # if query is contig_2
                elif query == gene2_contig.name:
                    link = CrossLink((contig_2_gene_content_track, query_start, query_end),
                                     (contig_1_gene_content_track, target_start, target_end),
                                     color=color,
                                     border=color,
                                     flip=False)
                    diagram.cross_track_links.append(link)

        ############################################### Draw and Export ################################################

        diagram.draw(format="linear",
                     orientation="landscape",
                     pagesize=(50 * cm, 25 * cm),
                     fragments=1,
                     start=0,
                     end=max_len)


        diagram.write('%s/%s.%s' % (path_to_output_act_folder, folder_name, flk_plot_fmt), flk_plot_fmt)


    # move plot to corresponding folder
    if match_category == 'end_match':
        os.system('mv %s/%s.%s %s/' % (path_to_output_act_folder, folder_name, flk_plot_fmt, pwd_at_ends_plot_folder))
    elif match_category == 'full_length_match':
        os.system('mv %s/%s.%s %s/' % (path_to_output_act_folder, folder_name, flk_plot_fmt, pwd_full_contig_match_plot_folder))
    else:
        os.system('mv %s/%s.%s %s/' % (path_to_output_act_folder, folder_name, flk_plot_fmt, pwd_normal_plot_folder))


    # remove temporary folder
    if keep_temp == 0:
        shutil.rmtree('%s/%s' % (path_to_output_act_folder, folder_name), ignore_errors=True)
        if os.path.isdir('%s/%s' % (path_to_output_act_folder, folder_name)):
            shutil.rmtree('%s/%s' % (path_to_output_act_folder, folder_name), ignore_errors=True)


def remove_bidirection(input_file, candidate2identity_dict, output_file):
    input = open(input_file)
    output = open(output_file, 'w')

    # get overall list
    overall = []
    for each in input:
        overall.append(each.strip())

    # get overlap list
    tmp_list = []
    overlap_list = []
    for each in overall:
        each_split = each.split('\t')
        each_reverse = '%s\t%s' % (each_split[1], each_split[0])
        tmp_list.append(each)
        if each_reverse in tmp_list:
            overlap_list.append(each)

    # get non-overlap list
    non_overlap_list = []
    for each in overall:
        each_split = each.split('\t')
        each_reverse = '%s\t%s' % (each_split[1], each_split[0])
        if (each not in overlap_list) and (each_reverse not in overlap_list):
            non_overlap_list.append(each)

    # get output
    for each in non_overlap_list:
        each_split = each.split('\t')
        each_concatenated = '%s___%s' % (each_split[0], each_split[1])
        output.write('%s\t%s\n' % (each, candidate2identity_dict[each_concatenated]))
    for each in overlap_list:
        each_concatenated = '%s___%s' % (each.split('\t')[0], each.split('\t')[1])
        output.write('%s\t%s\n' % (each, candidate2identity_dict[each_concatenated]))
    output.close()


def export_HGT_query_to_subjects(pwd_BM_HGTs, pwd_blast_subjects_in_one_line, pwd_query_to_subjects_file):

    HGT_candidates = set()
    for HGT_pair in open(pwd_BM_HGTs):
        HGT_pair_split = HGT_pair.strip().split('\t')
        gene_1 = HGT_pair_split[0]
        gene_2 = HGT_pair_split[1]
        HGT_candidates.add(gene_1)
        HGT_candidates.add(gene_2)

    query_subjects_dict = {}
    for each_gene in open(pwd_blast_subjects_in_one_line):
        each_gene_split = each_gene.strip().split('\t')
        query = each_gene_split[0].split('|')[1]
        subjects = [i.split('|')[1] for i in each_gene_split[1:]]
        if query in HGT_candidates:
            query_subjects_dict[query] = subjects

    pwd_query_to_subjects_file_handle = open(pwd_query_to_subjects_file, 'w')
    for each in query_subjects_dict:
        for_out = '%s\t%s\n' % (each, ','.join(query_subjects_dict[each]))
        pwd_query_to_subjects_file_handle.write(for_out)
    pwd_query_to_subjects_file_handle.close()


def filter_blast_results_worker(argument_list):
    pwd_blast_results = argument_list[0]
    align_len_cutoff = argument_list[1]
    cover_cutoff = argument_list[2]
    genome_name_list = argument_list[3]
    pwd_qualified_iden_file = argument_list[4]

    get_qualigied_blast_hits(pwd_blast_results, align_len_cutoff, cover_cutoff, genome_name_list, pwd_qualified_iden_file)


def get_g2g_identities_worker(argument_list):
    pwd_qualified_iden_file = argument_list[0]
    name_to_group_number_dict = argument_list[1]
    pwd_qualified_iden_file_g2g = argument_list[2]

    qualified_identities_g_g = open(pwd_qualified_iden_file_g2g, 'w')
    for each_identity in open(pwd_qualified_iden_file):
        each_identity_split = each_identity.strip().split('\t')
        query = each_identity_split[0]
        subject = each_identity_split[1]
        identity = float(each_identity_split[2])
        query_genome_name = '_'.join(query.split('_')[:-1])
        subject_genome_name = '_'.join(subject.split('_')[:-1])
        query_group = name_to_group_number_dict[query_genome_name].split('_')[0]
        subject_group = name_to_group_number_dict[subject_genome_name].split('_')[0]
        paired_group_list = [query_group, subject_group]

        # query and subjects name sorted by alphabet order here
        paired_group_list_sorted = sorted(paired_group_list)
        g_g = '%s_%s' % (paired_group_list_sorted[0], paired_group_list_sorted[1])
        qualified_identities_g_g.write('%s\t%s\n' % (g_g, identity))

    qualified_identities_g_g.close()


def get_HGT_worker(argument_list):
    pwd_qualified_iden_file = argument_list[0]
    name_to_group_number_dict = argument_list[1]
    pwd_qual_idens_with_group = argument_list[2]
    pwd_qual_idens_subjects_in_one_line = argument_list[3]
    pwd_hgt_candidates_with_group = argument_list[4]
    pwd_hgt_candidates_only_gene = argument_list[5]
    group_pair_iden_cutoff_dict = argument_list[6]


    file_path, file_basename, file_extension = sep_path_basename_ext(pwd_qual_idens_with_group)
    pwd_qual_idens_with_group_tmp = '%s/%s_tmp.%s' % (file_path, file_basename, file_extension)

    qualified_matches_with_group = open(pwd_qual_idens_with_group_tmp, 'w')
    for qualified_identity in open(pwd_qualified_iden_file):
        qualified_identity_split = qualified_identity.strip().split('\t')
        query = qualified_identity_split[0]
        query_split = query.split('_')
        query_bin = '_'.join(query_split[:-1])
        subject = qualified_identity_split[1]
        subject_split = subject.split('_')
        subject_bin = '_'.join(subject_split[:-1])
        identity = float(qualified_identity_split[2])
        file_write = '%s|%s\t%s|%s|%s\n' % (
            name_to_group_number_dict[query_bin], query, name_to_group_number_dict[subject_bin], subject, str(identity))
        qualified_matches_with_group.write(file_write)
    qualified_matches_with_group.close()

    # sort
    os.system('cat %s | sort > %s' % (pwd_qual_idens_with_group_tmp, pwd_qual_idens_with_group))
    os.system('rm %s' % pwd_qual_idens_with_group_tmp)

    # put subjects in one line
    get_hits_group(pwd_qual_idens_with_group, pwd_qual_idens_subjects_in_one_line)

    # get HGT candidates
    get_candidates(pwd_qual_idens_subjects_in_one_line,
                   pwd_hgt_candidates_with_group,
                   pwd_hgt_candidates_only_gene,
                   group_pair_iden_cutoff_dict)


def BM(args, config_dict):

    def do(plot_identity):
        current_group_pair_identities_array = np.array(current_group_pair_identities)
        current_group_pair_identity_cut_off = np.percentile(current_group_pair_identities_array, identity_percentile)
        current_group_pair_identity_cut_off = float("{0:.2f}".format(current_group_pair_identity_cut_off))
        current_group_pair_name_split = current_group_pair_name.split('_')
        current_group_pair_name_swapped = '%s_%s' % (current_group_pair_name_split[1], current_group_pair_name_split[0])
        group_pair_iden_cutoff_dict[current_group_pair_name] = current_group_pair_identity_cut_off
        group_pair_iden_cutoff_dict[current_group_pair_name_swapped] = current_group_pair_identity_cut_off

        if current_group_pair_name != current_group_pair_name_swapped:
            group_pair_iden_cutoff_file.write(
                '%s\t%s\n' % (current_group_pair_name, current_group_pair_identity_cut_off))

        # check length
        if len(current_group_pair_identities) >= minimum_plot_number:
            if current_group_pair_name == current_group_pair_name_swapped:
                if plot_identity is True:
                    plot_identity_list(current_group_pair_identities, 'None', current_group_pair_name, pwd_iden_distrib_plot_folder)
            else:
                if plot_identity is True:
                    plot_identity_list(current_group_pair_identities, current_group_pair_identity_cut_off, current_group_pair_name, pwd_iden_distrib_plot_folder)

            #report_and_log(("Plotting identity distribution (%dth): %s" % (ploted_group, current_group_pair_name)), pwd_log_file, keep_quiet)

        else:
            with open(pwd_unploted_groups_file, 'a') as unploted_groups_handle:
                unploted_groups_handle.write('%s\t%s\n' % (current_group_pair_name, len(current_group_pair_identities)))
            #report_and_log(("Plotting identity distribution (%dth): %s, blast hits < %d, skipped" % (ploted_group, current_group_pair_name, minimum_plot_number)), pwd_log_file, keep_quiet)

    output_prefix =             args['p']
    grouping_level =            args['r']
    grouping_file =             args['g']
    cover_cutoff =              args['cov']
    align_len_cutoff =          args['al']
    flanking_length_kbp =       args['flk']
    identity_percentile =       args['ip']
    end_match_identity_cutoff = args['ei']
    num_threads =               args['t']
    No_Eb_Check =               args['NoEbCheck']
    plot_identity =             args['plot_iden']
    keep_quiet =                args['quiet']
    keep_temp =                 args['tmp']


    # get path to current script
    pwd_blastn_exe = config_dict['blastn']

    flanking_length = flanking_length_kbp * 1000

    warnings.filterwarnings("ignore")

    #################################### find matched grouping file if not provided  ###################################

    if grouping_level is None:
        grouping_level = 'x'

    MetaCHIP_wd =   '%s_MetaCHIP_wd'              % output_prefix
    pwd_log_file =  '%s/%s_%s_BM_%s.log'  % (MetaCHIP_wd, output_prefix, grouping_level, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))


    pwd_grouping_file = ''
    group_num = 0
    if grouping_file is None:

        grouping_file_re = '%s/%s_%s*_grouping.txt' % (MetaCHIP_wd, output_prefix, grouping_level)
        grouping_file_list = [os.path.basename(file_name) for file_name in glob.glob(grouping_file_re)]

        if len(grouping_file_list) == 1:
            detected_grouping_file = grouping_file_list[0]
            pwd_grouping_file = '%s/%s' % (MetaCHIP_wd, detected_grouping_file)
            group_num = get_group_num_from_grouping_file(pwd_grouping_file)
            report_and_log(('Found grouping file %s, input genomes were clustered into %s groups' % (detected_grouping_file, group_num)), pwd_log_file, keep_quiet)

        elif len(grouping_file_list) == 0:
            report_and_log(('No grouping file detected, please specify with "-g" option'), pwd_log_file, keep_quiet)
            exit()

        else:
            report_and_log(('Multiple grouping file detected, please specify with "-g" option'), pwd_log_file, keep_quiet)
            exit()

    else:  # with provided grouping file
        pwd_grouping_file = grouping_file
        group_num = get_group_num_from_grouping_file(pwd_grouping_file)


    ############################################# define file/folder names #############################################

    MetaCHIP_op_folder = '%s_%s%s_HGTs_ip%s_al%sbp_c%s_ei%sbp_f%skbp' % (output_prefix, grouping_level, group_num, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)

    blast_result_folder =                               '%s_all_blastn_results'                           % (output_prefix)
    combined_ffn_file =                                 '%s_all_combined_ffn.fasta'                       % (output_prefix)
    blast_result_filtered_folder =                      '%s_%s%s_blastn_results_filtered'                 % (output_prefix, grouping_level, group_num)
    blast_result_filtered_folder_g2g =                  '%s_%s%s_1_blastn_results_filtered_g2g'           % (output_prefix, grouping_level, group_num)
    blast_result_filtered_folder_with_group =           '%s_%s%s_2_blastn_results_filtered_with_group'    % (output_prefix, grouping_level, group_num)
    blast_result_filtered_folder_in_one_line =          '%s_%s%s_3_blastn_results_filtered_in_one_line'   % (output_prefix, grouping_level, group_num)
    op_candidates_only_gene_folder_name =               '%s_%s%s_4_HGTs_only_id'                          % (output_prefix, grouping_level, group_num)
    op_candidates_with_group_folder_name =              '%s_%s%s_5_HGTs_with_group'                       % (output_prefix, grouping_level, group_num)
    gbk_folder =                                        '%s_%s%s_gbk_files'                               % (output_prefix, grouping_level, group_num)
    iden_distrib_plot_folder =                          '%s_%s%s_identity_distribution'                   % (output_prefix, grouping_level, group_num)
    qual_idens_file =                                   '%s_%s%s_blastn_results_filtered.tab'             % (output_prefix, grouping_level, group_num)
    qual_idens_file_gg =                                '%s_%s%s_qualified_iden_gg.txt'                   % (output_prefix, grouping_level, group_num)
    qual_idens_file_gg_sorted =                         '%s_%s%s_qualified_iden_gg_sorted.txt'            % (output_prefix, grouping_level, group_num)
    subjects_in_one_line_filename =                     '%s_%s%s_subjects_in_one_line.txt'                % (output_prefix, grouping_level, group_num)
    HGT_query_to_subjects_filename =                    '%s_%s%s_HGT_query_to_subjects.txt'               % (output_prefix, grouping_level, group_num)
    group_pair_iden_cutoff_file_name =                  '%s_%s%s_identity_cutoff.txt'                     % (output_prefix, grouping_level, group_num)
    op_candidates_with_group_file_name =                '%s_%s%s_HGTs_with_group.txt'                     % (output_prefix, grouping_level, group_num)
    op_candidates_only_gene_file_name =                 '%s_%s%s_HGTs_only_id.txt'                        % (output_prefix, grouping_level, group_num)
    op_candidates_only_gene_file_name_uniq =            '%s_%s%s_HGTs_uniq.txt'                           % (output_prefix, grouping_level, group_num)
    op_candidates_BM =                                  '%s_%s%s_HGTs_BM.txt'                             % (output_prefix, grouping_level, group_num)
    op_candidates_seq_nc =                              '%s_%s%s_HGTs_BM_nc.fasta'                        % (output_prefix, grouping_level, group_num)
    op_act_folder_name =                                '%s_%s%s_Flanking_region_plots'                   % (output_prefix, grouping_level, group_num)
    grouping_file_with_id_filename =                    '%s_%s%s_grouping_with_id.txt'                    % (output_prefix, grouping_level, group_num)
    unploted_groups_file =                              '0_unploted_groups.txt'
    normal_folder_name =                                '1_Plots_normal'
    end_match_folder_name =                             '2_Plots_end_match'
    full_length_match_folder_name =                     '3_Plots_full_length_match'


    pwd_MetaCHIP_op_folder =                       '%s/%s'       % (MetaCHIP_wd, MetaCHIP_op_folder)
    pwd_gbk_folder =                               '%s/%s'       % (MetaCHIP_wd, gbk_folder)
    pwd_combined_ffn_file =                        '%s/%s'       % (MetaCHIP_wd, combined_ffn_file)
    pwd_blast_result_folder =                      '%s/%s'       % (MetaCHIP_wd, blast_result_folder)
    pwd_blast_result_filtered_folder =             '%s/%s'       % (MetaCHIP_wd, blast_result_filtered_folder)
    pwd_qualified_iden_file =                      '%s/%s'       % (MetaCHIP_wd, qual_idens_file)
    pwd_blast_result_filtered_folder_g2g =         '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, blast_result_filtered_folder_g2g)
    pwd_blast_result_filtered_folder_with_group =  '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, blast_result_filtered_folder_with_group)
    pwd_blast_result_filtered_folder_in_one_line = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, blast_result_filtered_folder_in_one_line)
    pwd_iden_distrib_plot_folder =                 '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, iden_distrib_plot_folder)
    pwd_qual_iden_file_gg =                        '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, qual_idens_file_gg)
    pwd_qual_iden_file_gg_sorted =                 '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, qual_idens_file_gg_sorted)
    pwd_unploted_groups_file =                     '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, iden_distrib_plot_folder, unploted_groups_file)
    pwd_subjects_in_one_line =                     '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, subjects_in_one_line_filename)
    pwd_HGT_query_to_subjects_file =               '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, HGT_query_to_subjects_filename)
    pwd_group_pair_iden_cutoff_file =              '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, group_pair_iden_cutoff_file_name)
    pwd_op_candidates_with_group_file =            '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_with_group_file_name)
    pwd_op_candidates_only_gene_file =             '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_only_gene_file_name)
    pwd_op_candidates_with_group_folder =          '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_with_group_folder_name)
    pwd_op_candidates_only_gene_folder =           '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_only_gene_folder_name)
    pwd_op_candidates_only_gene_file_uniq =        '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_only_gene_file_name_uniq)
    pwd_op_candidates_BM =                         '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_BM)
    pwd_op_candidates_seq_nc =                     '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_seq_nc)
    pwd_op_act_folder =                            '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name)
    pwd_grouping_file_with_id =                    '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, grouping_file_with_id_filename)
    pwd_normal_folder =                            '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name, normal_folder_name)
    pwd_end_match_folder =                         '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name, end_match_folder_name)
    pwd_full_length_match_folder =                 '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name, full_length_match_folder_name)


    ####################################################################################################################

    # create outputs folder
    force_create_folder(pwd_MetaCHIP_op_folder)

    # index grouping file
    index_grouping_file(pwd_grouping_file, pwd_grouping_file_with_id)

    # create genome_group_dict and genome_list
    genome_name_list = []
    name_to_group_number_dict = {}
    name_to_group_dict = {}
    for each_bin in open(pwd_grouping_file_with_id):
        each_bin_split = each_bin.strip().split(',')
        bin_name = each_bin_split[1]
        bin_group_number = each_bin_split[0]
        bin_group = bin_group_number.split('_')[0]
        genome_name_list.append(bin_name)
        name_to_group_number_dict[bin_name] = bin_group_number
        name_to_group_dict[bin_name] = bin_group


    ############################################### filter blastn results ##############################################

    # check whether previous filter results exist
    filter_blastn_results = 0
    if os.path.isdir(pwd_blast_result_filtered_folder) is True:
        blast_result_filtered_file_re = '%s/*_filtered.tab' % pwd_blast_result_filtered_folder
        blast_result_filtered_file_list = [os.path.basename(file_name) for file_name in glob.glob(blast_result_filtered_file_re)]
        if len(blast_result_filtered_file_list) != len(genome_name_list):
            report_and_log(('A difference in numbers between the filtered blast results and qualified genomes were found, will re-run the blast results filter step for all qualified genomes'), pwd_log_file, keep_quiet)
            filter_blastn_results = 1
    else:
        filter_blastn_results = 1

    # filter blast results if filtered blast results not exist
    if filter_blastn_results == 0:
        report_and_log(('Filtered blastn results at specified taxonomic rank detected from folder %s. HGT analysis will be performed based on these files.' % blast_result_filtered_folder), pwd_log_file, keep_quiet)

    else:
        report_and_log(('Filtering blast matches with the following criteria: Query genome != Subject genome, Alignment length >= %sbp and coverage >= %s%s' % (align_len_cutoff, cover_cutoff, '%')), pwd_log_file, keep_quiet)

        # create folder
        force_create_folder(pwd_blast_result_filtered_folder)

        # filter blastn results with multiprocessing
        blast_result_file_re = '%s/*_blastn.tab' % pwd_blast_result_folder
        blast_result_file_list = [os.path.basename(file_name) for file_name in glob.glob(blast_result_file_re)]
        if len(blast_result_file_list) == 0:
            report_and_log(('No blast results detected, program exited!'), pwd_log_file, keep_quiet)
            exit()

        list_for_multiple_arguments_filter_blast_results = []
        for blast_result_file in blast_result_file_list:
            genome_name = blast_result_file.split('_blastn')[0]
            if genome_name in genome_name_list:
                pwd_blast_result_file = '%s/%s' % (pwd_blast_result_folder, blast_result_file)
                blast_result_filtered_file = '%s_filtered.tab' % ('.'.join(blast_result_file.split('.')[:-1]))
                pwd_blast_result_filtered_file = '%s/%s' % (pwd_blast_result_filtered_folder, blast_result_filtered_file)
                list_for_multiple_arguments_filter_blast_results.append([pwd_blast_result_file, align_len_cutoff, cover_cutoff, genome_name_list, pwd_blast_result_filtered_file])

        # filter_blast_results with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(filter_blast_results_worker, list_for_multiple_arguments_filter_blast_results)
        pool.close()
        pool.join()


    # combine qualified blast hits
    report_and_log(('Combining filtered blastn results'), pwd_log_file, keep_quiet)
    os.system('cat %s/*_blastn_filtered.tab > %s' % (pwd_blast_result_filtered_folder, pwd_qualified_iden_file))


    ########################################### get group-to-group identities ##########################################

    report_and_log(('Get group-to-group identities with %s cores' % num_threads), pwd_log_file, keep_quiet)

    # create folder
    force_create_folder(pwd_blast_result_filtered_folder_g2g)

    # get file list
    blast_result_filtered_file_re = '%s/*_filtered.tab' % pwd_blast_result_filtered_folder
    blast_result_filtered_file_list = [os.path.basename(file_name) for file_name in glob.glob(blast_result_filtered_file_re)]

    list_for_multiple_arguments_get_g2g_identities = []
    for filtered_blast_result in blast_result_filtered_file_list:
        pwd_filtered_blast_result = '%s/%s' % (pwd_blast_result_filtered_folder, filtered_blast_result)
        pwd_filtered_blast_result_g2g = '%s/%s_g2g.tab' % (pwd_blast_result_filtered_folder_g2g, filtered_blast_result.split('.tab')[0])
        list_for_multiple_arguments_get_g2g_identities.append([pwd_filtered_blast_result, name_to_group_number_dict, pwd_filtered_blast_result_g2g])

    # get group-to-group identities with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(get_g2g_identities_worker, list_for_multiple_arguments_get_g2g_identities)
    pool.close()
    pool.join()

    # combine g2g_files and sort it
    os.system('cat %s/*_g2g.tab > %s' % (pwd_blast_result_filtered_folder_g2g, pwd_qual_iden_file_gg))
    os.system('cat %s | sort > %s' % (pwd_qual_iden_file_gg, pwd_qual_iden_file_gg_sorted))


    ############ plot identity distribution between groups and get cutoff according to specified percentile ############

    # create folder to hold group-group identity distribution plot
    report_and_log(('Plotting identity distribution between each pair of groups'), pwd_log_file, keep_quiet)
    os.makedirs(pwd_iden_distrib_plot_folder)

    # get identities for each group pair, plot identity distribution and generate group_pair to identity dict
    with open(pwd_unploted_groups_file, 'a') as unploted_groups_handle:
        unploted_groups_handle.write('Group\tHits_number\n')
    current_group_pair_name = ''
    current_group_pair_identities = []
    group_pair_iden_cutoff_dict = {}
    ploted_group = 1
    minimum_plot_number = 10
    group_pair_iden_cutoff_file = open(pwd_group_pair_iden_cutoff_file, 'w')
    for each_identity_g in open(pwd_qual_iden_file_gg_sorted):
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
                do(plot_identity)
                ploted_group += 1
                # restore current_group_pair_name and current_group_pair_identities for next group pair
                current_group_pair_name = group_pair
                current_group_pair_identities = []
                current_group_pair_identities.append(identity)

    # for the last group
    do(plot_identity)
    group_pair_iden_cutoff_file.close()


    ############################### add group to blast hits and put subjects in one line ###############################

    report_and_log(('Analyzing Blast hits to get HGT candidates with %s cores' % num_threads), pwd_log_file, keep_quiet)

    # create folder
    force_create_folder(pwd_blast_result_filtered_folder_with_group)
    force_create_folder(pwd_blast_result_filtered_folder_in_one_line)
    force_create_folder(pwd_op_candidates_with_group_folder)
    force_create_folder(pwd_op_candidates_only_gene_folder)

    list_for_multiple_arguments_get_HGT = []
    for filtered_blast_result in blast_result_filtered_file_list:
        genome_id = filtered_blast_result.split('_blastn_filtered.tab')[0]
        pwd_filtered_blast_result = '%s/%s' % (pwd_blast_result_filtered_folder, filtered_blast_result)
        pwd_filtered_blast_result_with_group =  '%s/%s_with_group_sorted.tab'    % (pwd_blast_result_filtered_folder_with_group, filtered_blast_result.split('.tab')[0])
        pwd_filtered_blast_result_in_one_line = '%s/%s_subjects_in_one_line.tab' % (pwd_blast_result_filtered_folder_in_one_line, filtered_blast_result.split('.tab')[0])
        pwd_hgt_candidates_with_group =         '%s/%s_HGTs_with_group.txt'      % (pwd_op_candidates_with_group_folder, genome_id)
        pwd_hgt_candidates_only_gene =          '%s/%s_HGTs_only_gene.txt'       % (pwd_op_candidates_only_gene_folder, genome_id)

        list_for_multiple_arguments_get_HGT.append([pwd_filtered_blast_result,
                                                    name_to_group_number_dict,
                                                    pwd_filtered_blast_result_with_group,
                                                    pwd_filtered_blast_result_in_one_line,
                                                    pwd_hgt_candidates_with_group,
                                                    pwd_hgt_candidates_only_gene,
                                                    group_pair_iden_cutoff_dict])

    # add group to blast hits with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(get_HGT_worker, list_for_multiple_arguments_get_HGT)
    pool.close()
    pool.join()


    ################################ remove bidirection and add identity to output file ################################

    # combine pwd_op_candidates_with_group_file and pwd_op_candidates_only_gene_file
    os.system('cat %s/*_HGTs_with_group.txt > %s' % (pwd_op_candidates_with_group_folder, pwd_op_candidates_with_group_file))
    os.system('cat %s/*_HGTs_only_gene.txt > %s'  % (pwd_op_candidates_only_gene_folder, pwd_op_candidates_only_gene_file))


    # remove bidirection and add identity to output file
    candidate2identity_dict = {}
    for each_candidate_pair in open(pwd_op_candidates_with_group_file):
        each_candidate_pair_split = each_candidate_pair.strip().split('\t')
        recipient_gene = each_candidate_pair_split[0].split('|')[1]
        donor_gene = each_candidate_pair_split[1].split('|')[1]
        identity = float(each_candidate_pair_split[1].split('|')[2])
        candidate2identity_key = '%s___%s' % (recipient_gene, donor_gene)
        candidate2identity_dict[candidate2identity_key] = identity

    remove_bidirection(pwd_op_candidates_only_gene_file, candidate2identity_dict, pwd_op_candidates_only_gene_file_uniq)


    ############################################### plot flanking region ###############################################

    report_and_log(('Plotting flanking regions with %s cores' % num_threads), pwd_log_file, keep_quiet)

    # create folder to hold ACT output
    os.makedirs(pwd_op_act_folder)
    os.makedirs(pwd_normal_folder)
    os.makedirs(pwd_end_match_folder)
    os.makedirs(pwd_full_length_match_folder)

    # initialize manager.dict
    manager = mp.Manager()
    candidates_2_contig_match_category_dict_mp = manager.dict()

    list_for_multiple_arguments_flanking_regions = []
    for match in open(pwd_op_candidates_only_gene_file_uniq):
        list_for_multiple_arguments_flanking_regions.append([match, pwd_gbk_folder, flanking_length, align_len_cutoff, name_to_group_number_dict, pwd_op_act_folder,
                                                             pwd_normal_folder, pwd_end_match_folder, pwd_full_length_match_folder, pwd_blastn_exe, keep_temp,
                                                             candidates_2_contig_match_category_dict_mp, end_match_identity_cutoff, No_Eb_Check])

    pool_flanking_regions = mp.Pool(processes=num_threads)
    pool_flanking_regions.map(get_gbk_blast_act2, list_for_multiple_arguments_flanking_regions)
    pool_flanking_regions.close()
    pool_flanking_regions.join()

    # convert mp.dict to normal dict
    candidates_2_contig_match_category_dict = {each_key: each_value for each_key, each_value in candidates_2_contig_match_category_dict_mp.items()}


    ################################################ get BM output file ################################################

    # add at_end information to output file
    BM_output_file_handle = open(pwd_op_candidates_BM, 'w')
    BM_output_file_handle.write('Gene_1\tGene_2\tGene_1_group\tGene_2_group\tIdentity\tend_match\tfull_length_match\n')
    for each_candidate in open(pwd_op_candidates_only_gene_file_uniq):
        each_candidate_split = each_candidate.strip().split('\t')
        recipient_gene = each_candidate_split[0]
        recipient_genome = '_'.join(recipient_gene.split('_')[:-1])
        recipient_genome_group_id = name_to_group_number_dict[recipient_genome]
        recipient_genome_group = recipient_genome_group_id.split('_')[0]
        donor_gene = each_candidate_split[1]
        donor_genome = '_'.join(donor_gene.split('_')[:-1])
        donor_genome_group_id = name_to_group_number_dict[donor_genome]
        donor_genome_group = donor_genome_group_id.split('_')[0]
        identity = each_candidate_split[2]
        concatenated = '%s___%s' % (recipient_gene, donor_gene)

        end_match_value = 'no'
        if candidates_2_contig_match_category_dict[concatenated] == 'end_match':
            end_match_value = 'yes'

        full_length_match_value = 'no'
        if candidates_2_contig_match_category_dict[concatenated] == 'full_length_match':
            full_length_match_value = 'yes'

        # write to output files
        BM_output_file_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_group, donor_genome_group, identity, end_match_value, full_length_match_value))
    BM_output_file_handle.close()


    ####################################### export gene clusters for PG approach #######################################

    os.system('cat %s/*_in_one_line.tab > %s' % (pwd_blast_result_filtered_folder_in_one_line, pwd_subjects_in_one_line))

    export_HGT_query_to_subjects(pwd_op_candidates_BM, pwd_subjects_in_one_line, pwd_HGT_query_to_subjects_file)


    ################################### export nc and aa sequence of predicted HGTs ####################################

    report_and_log(('Extracting nc sequences for BM predicted HGTs'), pwd_log_file, keep_quiet)

    # get qualified HGT candidates
    HGT_candidates_qualified = set()
    for each_candidate_2 in open(pwd_op_candidates_BM):
        if not each_candidate_2.startswith('Gene_1\tGene_2\tGene_1_group'):
            each_candidate_2_split = each_candidate_2.strip().split('\t')
            end_match = each_candidate_2_split[5]
            full_length_match = each_candidate_2_split[6]
            if (end_match == 'no') and (full_length_match == 'no'):
                HGT_candidates_qualified.add(each_candidate_2_split[0])
                HGT_candidates_qualified.add(each_candidate_2_split[1])

    candidates_seq_nc_handle = open(pwd_op_candidates_seq_nc, 'w')
    for each_seq in SeqIO.parse(pwd_combined_ffn_file, 'fasta'):
        if each_seq.id in HGT_candidates_qualified:
            SeqIO.write(each_seq, candidates_seq_nc_handle, 'fasta')
    candidates_seq_nc_handle.close()


    ############################################## remove temporary files ##############################################

    if keep_temp == 0:
        report_and_log(('Deleting temporary files'), pwd_log_file, keep_quiet)
        os.remove(pwd_qualified_iden_file)
        os.remove(pwd_qual_iden_file_gg)
        os.remove(pwd_qual_iden_file_gg_sorted)
        os.remove(pwd_subjects_in_one_line)
        os.remove(pwd_op_candidates_only_gene_file_uniq)
        os.remove(pwd_op_candidates_with_group_file)
        os.remove(pwd_op_candidates_only_gene_file)
        # os.remove(pwd_HGT_query_to_subjects_file) need this file in the PG approach
        os.system('rm -r %s' % pwd_blast_result_filtered_folder_g2g)
        os.system('rm -r %s' % pwd_blast_result_filtered_folder_with_group)
        os.system('rm -r %s' % pwd_blast_result_filtered_folder_in_one_line)
        os.system('rm -r %s' % pwd_op_candidates_only_gene_folder)
        os.system('rm -r %s' % pwd_op_candidates_with_group_folder)

    # report
    report_and_log(('Done for Best-match approach!'), pwd_log_file, keep_quiet)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',             required=True,  help='output prefix')
    parser.add_argument('-r',             required=False, default=None, help='grouping rank')
    parser.add_argument('-g',             required=False, default=None, help='grouping file')
    parser.add_argument('-cov',           required=False, type=int, default=75, help='coverage cutoff, default: 75')
    parser.add_argument('-al',            required=False, type=int, default=200, help='alignment length cutoff, default: 200')
    parser.add_argument('-flk',           required=False, type=int, default=10, help='the length of flanking sequences to plot (Kbp), default: 10')
    parser.add_argument('-ip',            required=False, type=int, default=90, help='identity percentile cutoff, default: 90')
    parser.add_argument('-ei',            required=False, type=float, default=90, help='end match identity cutoff, default: 95')
    parser.add_argument('-t',             required=False, type=int, default=1, help='number of threads, default: 1')
    parser.add_argument('-plot_iden',     required=False, action="store_true", help='plot identity distribution')
    parser.add_argument('-NoEbCheck',     required=False, action="store_true", help='disable end break and contig match check for fast processing, not recommend for metagenome-assembled genomes (MAGs)')
    parser.add_argument('-force',         required=False, action="store_true", help='overwrite previous results')
    parser.add_argument('-quiet',         required=False, action="store_true", help='Do not report progress')
    parser.add_argument('-tmp',           required=False, action="store_true", help='keep temporary files')

    args = vars(parser.parse_args())


    BM(args, config_dict)

