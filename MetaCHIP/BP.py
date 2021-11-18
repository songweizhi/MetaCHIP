#!/usr/bin/env python3

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
import re
import sys
import copy
import glob
import shutil
import platform
import warnings
import argparse
import itertools
import subprocess
import numpy as np
import multiprocessing as mp
from time import sleep
from ete3 import Tree
from Bio import SeqIO, AlignIO, Align, Phylo
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm
from datetime import datetime
from string import ascii_uppercase
from scipy.stats import gaussian_kde
from distutils.spawn import find_executable
from MetaCHIP.MetaCHIP_config import config_dict


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


def force_create_folder(folder_to_create):

    rm_rd = 0
    while os.path.isdir(folder_to_create) is True:
        shutil.rmtree(folder_to_create, ignore_errors=True)

        if rm_rd >= 10:
            print('Failed in removing %s, program exited!' % folder_to_create)
            exit()

        rm_rd += 1
        sleep(1)

    os.mkdir(folder_to_create)


def rm_folder_file(target_re):
    target_list = glob.glob(target_re)

    for target in target_list:

        if os.path.isdir(target) is True:
            os.system('rm -r %s' % target)

        elif os.path.isfile(target) is True:
            os.system('rm %s' % target)


def filter_blast_results_worker(argument_list):
    pwd_blast_results = argument_list[0]
    align_len_cutoff = argument_list[1]
    cover_cutoff = argument_list[2]
    pwd_qualified_iden_file = argument_list[3]

    pwd_qualified_iden_file_handle = open(pwd_qualified_iden_file, 'w')
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

        if (align_len >= int(align_len_cutoff)) and (query_bin_name != subject_bin_name) and (coverage_q >= cover_cutoff) and (coverage_s >= cover_cutoff):
            pwd_qualified_iden_file_handle.write(match)

    pwd_qualified_iden_file_handle.close()

    # remove if empty
    if os.stat(pwd_qualified_iden_file).st_size == 0:
        os.system('rm %s' % pwd_qualified_iden_file)


def get_number_of_group(grouping_file):

    group_list = []
    for each_genome in open(grouping_file):
        each_genome_split = each_genome.strip().split(',')
        group_id = each_genome_split[0]
        if group_id not in group_list:
            group_list.append(group_id)

    return len(group_list)


def get_g2g_identities_worker(argument_list):
    pwd_qualified_iden_file     = argument_list[0]
    qualified_genome_list            = argument_list[1]
    name_to_group_number_dict   = argument_list[2]
    pwd_qualified_iden_file_g2g = argument_list[3]

    qualified_identities_g_g = open(pwd_qualified_iden_file_g2g, 'w')
    for each_identity in open(pwd_qualified_iden_file):
        each_identity_split = each_identity.strip().split('\t')
        query = each_identity_split[0]
        subject = each_identity_split[1]
        identity = float(each_identity_split[2])
        query_genome_name = '_'.join(query.split('_')[:-1])
        subject_genome_name = '_'.join(subject.split('_')[:-1])

        if (query_genome_name in qualified_genome_list) and (subject_genome_name in qualified_genome_list):
            query_group = name_to_group_number_dict[query_genome_name].split('_')[0]
            subject_group = name_to_group_number_dict[subject_genome_name].split('_')[0]
            paired_group_list = [query_group, subject_group]

            # query and subjects name sorted by alphabet order here
            paired_group_list_sorted = sorted(paired_group_list)
            g_g = '%s_%s' % (paired_group_list_sorted[0], paired_group_list_sorted[1])
            qualified_identities_g_g.write('%s\t%s\n' % (g_g, identity))

    qualified_identities_g_g.close()


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


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


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

        # if only one non-self subject was found, no matter which group it comes from, ignored. proceed if more than 1 non-self subject was found:
        if len(subjects_list) > 1:
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
                    sg_average = sg_sum/float(sg_subject_number)

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
                    nsg_average = nsg_sum/float(nsg_subject_number)

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
                    sg_average = sg_sum/float(sg_subject_number)

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
                            each_nsg_subject_iden = float(each_nsg_subject.split('|')[2])
                            each_nsg_subject_group = each_nsg_subject.split('|')[0].split('_')[0]
                            if each_nsg_subject_group == each_nsg:
                                if each_nsg_subject_iden > nsg_maximum:
                                    nsg_maximum = each_nsg_subject_iden
                                    nsg_maximum_gene = each_nsg_subject
                                nsg_sum += each_nsg_subject_iden
                                nsg_subject_number += 1
                        nsg_average = nsg_sum / float(nsg_subject_number)
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


def get_HGT_worker(argument_list):
    pwd_qualified_iden_file = argument_list[0]
    name_to_group_number_dict = argument_list[1]
    pwd_qual_idens_with_group = argument_list[2]
    pwd_qual_idens_subjects_in_one_line = argument_list[3]
    pwd_hgt_candidates_with_group = argument_list[4]
    pwd_hgt_candidates_only_gene = argument_list[5]
    group_pair_iden_cutoff_dict = argument_list[6]
    qualified_genome_list = argument_list[7]

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
        if (query_bin in qualified_genome_list) and (subject_bin in qualified_genome_list):
            file_write = '%s|%s\t%s|%s|%s\n' % (name_to_group_number_dict[query_bin], query, name_to_group_number_dict[subject_bin], subject, str(identity))
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


class BinRecord(object):

    def __init__(self, name, group, group_without_underscore):
        self.name = name
        self.group = group
        self.group_without_underscore = group_without_underscore


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
    query_cov_total = query_matched_len_total/float(query_len)
    subject_cov_total = subject_matched_len_total/float(subject_len)

    # get match category
    match_category = 'normal'
    best_hit_end_gap_len = 200
    gap_cutoff_for_concatenating = 300

    # full length match: coverage cutoff 90%
    if (query_cov_total >= 0.9) or (subject_cov_total >= 0.9):
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
                            pass

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
            feature.qualifiers['locus_tag'][0] = 'Gene_%s' % bin_name_gbk_split[-1]

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
            feature_set.add_feature(feature, color=color, label=True, sigil='ARROW', arrowshaft_height=0.5, arrowhead_length=0.4, label_color=label_color, label_size=label_size, label_angle=label_angle, label_position="middle")


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


    ############################## check whether full length or end match ##############################

    # get match category
    if No_Eb_Check is True:
        match_category = 'normal'
        candidates_2_contig_match_category_dict[folder_name] = match_category
    else:

        # run blast
        os.system(command_blast_full_len)

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
            candidates_2_contig_match_category_dict[folder_name] = match_category
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


def subset_tree(tree_file_in, leaf_node_list, tree_file_out):
    tree_in = Tree(tree_file_in, format=0)
    tree_in.prune(leaf_node_list, preserve_branch_length=True)
    tree_in.write(format=0, outfile=tree_file_out)


def extract_gene_tree_seq_worker(argument_list):

    each_to_process =               argument_list[0]
    pwd_tree_folder =               argument_list[1]
    pwd_combined_faa_file_subset =  argument_list[2]
    pwd_blastp_exe =                argument_list[3]
    pwd_mafft_exe =                 argument_list[4]
    pwd_fasttree_exe =              argument_list[5]
    genome_to_group_dict =          argument_list[6]
    genome_name_list =              argument_list[7]
    HGT_query_to_subjects_dict =    argument_list[8]
    pwd_SCG_tree_all =              argument_list[9]

    gene_1 = each_to_process[0]
    gene_2 = each_to_process[1]
    HGT_genome_1 = '_'.join(gene_1.split('_')[:-1])
    HGT_genome_2 = '_'.join(gene_2.split('_')[:-1])
    paired_groups = [genome_to_group_dict[HGT_genome_1], genome_to_group_dict[HGT_genome_2]]


    each_to_process_concate = '___'.join(each_to_process)
    blast_output =            '%s/%s___%s_gene_tree_blast.tab'        % (pwd_tree_folder, gene_1, gene_2)
    blast_output_sorted =     '%s/%s___%s_gene_tree_blast_sorted.tab' % (pwd_tree_folder, gene_1, gene_2)
    gene_tree_seq =           '%s/%s___%s_gene_tree.seq'              % (pwd_tree_folder, gene_1, gene_2)
    gene_tree_seq_uniq =      '%s/%s___%s_gene_tree_uniq.seq'         % (pwd_tree_folder, gene_1, gene_2)
    self_seq =                '%s/%s___%s_gene_tree_selfseq.seq'      % (pwd_tree_folder, gene_1, gene_2)
    non_self_seq =            '%s/%s___%s_gene_tree_nonselfseq.seq'   % (pwd_tree_folder, gene_1, gene_2)
    pwd_seq_file_1st_aln =    '%s/%s___%s_gene_tree.1.aln'            % (pwd_tree_folder, gene_1, gene_2)
    pwd_seq_file_2nd_aln =    '%s/%s___%s_gene_tree.2.aln'            % (pwd_tree_folder, gene_1, gene_2)
    pwd_gene_tree_newick =    '%s/%s___%s_gene_tree.newick'           % (pwd_tree_folder, gene_1, gene_2)
    pwd_species_tree_newick = '%s/%s___%s_species_tree.newick'        % (pwd_tree_folder, gene_1, gene_2)

    ################################################## Get gene tree ###################################################

    current_gene_member_BM = set()
    current_gene_member_BM.add(gene_1)
    current_gene_member_BM.add(gene_2)

    if gene_1 in HGT_query_to_subjects_dict:
        for gene_1_subject in HGT_query_to_subjects_dict[gene_1]:
            current_gene_member_BM.add(gene_1_subject)
    if gene_2 in HGT_query_to_subjects_dict:
        for gene_2_subject in HGT_query_to_subjects_dict[gene_2]:
            current_gene_member_BM.add(gene_2_subject)

    current_gene_member_grouped = []
    for gene_member in current_gene_member_BM:
        gene_member_genome = '_'.join(gene_member.split('_')[:-1])
        if gene_member_genome in genome_name_list:
            current_gene_member_grouped.append(gene_member)

    current_gene_member_grouped_from_paired_group = []
    for gene_member in current_gene_member_grouped:
        current_gene_genome = '_'.join(gene_member.split('_')[:-1])
        current_genome_group = genome_to_group_dict[current_gene_genome]
        if current_genome_group in paired_groups:
            current_gene_member_grouped_from_paired_group.append(gene_member)

    # genes to extract
    if len(current_gene_member_grouped_from_paired_group) < 3:
        genes_to_extract_list = current_gene_member_grouped
    else:
        genes_to_extract_list = current_gene_member_grouped_from_paired_group

    # get sequences of othorlog group to build gene tree
    output_handle = open(gene_tree_seq, "w")
    extracted_gene_set = set()
    for seq_record in SeqIO.parse(pwd_combined_faa_file_subset, 'fasta'):
        # if seq_record.id in current_gene_member:
        if seq_record.id in genes_to_extract_list:
            output_handle.write('>%s\n' % seq_record.id)
            output_handle.write('%s\n' % str(seq_record.seq))
            extracted_gene_set.add(seq_record.id)
    output_handle.close()

    if (gene_1 in extracted_gene_set) and (gene_2 in extracted_gene_set):
        self_seq_handle = open(self_seq, 'w')
        non_self_seq_handle = open(non_self_seq, 'w')
        non_self_seq_num = 0
        for each_seq in SeqIO.parse(gene_tree_seq, 'fasta'):
            each_seq_genome_id = '_'.join(each_seq.id.split('_')[:-1])
            if each_seq.id in each_to_process:
                SeqIO.write(each_seq, self_seq_handle, 'fasta')
            elif each_seq_genome_id not in [HGT_genome_1, HGT_genome_2]:
                SeqIO.write(each_seq, non_self_seq_handle, 'fasta')
                non_self_seq_num += 1
        self_seq_handle.close()
        non_self_seq_handle.close()


        # run blast
        genome_subset = set()
        if non_self_seq_num > 0:
            os.system('%s -query %s -subject %s -outfmt 6 -out %s' % (pwd_blastp_exe, self_seq, non_self_seq, blast_output))
            os.system('cat %s | sort > %s' % (blast_output, blast_output_sorted))

            # get best match from each genome
            current_query_subject_genome = ''
            current_bit_score = 0
            current_best_match = ''
            best_match_list = []
            for each_hit in open(blast_output_sorted):
                each_hit_split = each_hit.strip().split('\t')
                query = each_hit_split[0]
                subject = each_hit_split[1]
                subject_genome = '_'.join(subject.split('_')[:-1])
                query_subject_genome = '%s___%s' % (query, subject_genome)
                bit_score = float(each_hit_split[11])
                if current_query_subject_genome == '':
                    current_query_subject_genome = query_subject_genome
                    current_bit_score = bit_score
                    current_best_match = subject
                elif current_query_subject_genome == query_subject_genome:
                    if bit_score > current_bit_score:
                        current_bit_score = bit_score
                        current_best_match = subject
                elif current_query_subject_genome != query_subject_genome:
                    best_match_list.append(current_best_match)
                    current_query_subject_genome = query_subject_genome
                    current_bit_score = bit_score
                    current_best_match = subject
            best_match_list.append(current_best_match)

            # export sequences
            gene_tree_seq_all = best_match_list + each_to_process
            gene_tree_seq_uniq_handle = open(gene_tree_seq_uniq, 'w')
            for each_seq2 in SeqIO.parse(gene_tree_seq, 'fasta'):
                if each_seq2.id in gene_tree_seq_all:
                    gene_tree_seq_uniq_handle.write('>%s\n' % each_seq2.id)
                    gene_tree_seq_uniq_handle.write('%s\n' % str(each_seq2.seq))
            gene_tree_seq_uniq_handle.close()

            cmd_mafft = '%s --quiet %s > %s' % (pwd_mafft_exe, gene_tree_seq_uniq, pwd_seq_file_1st_aln)
            for each_gene in SeqIO.parse(gene_tree_seq_uniq, 'fasta'):
                each_gene_genome = '_'.join(str(each_gene.id).split('_')[:-1])
                genome_subset.add(each_gene_genome)
        else:
            cmd_mafft = '%s --quiet %s > %s' % (pwd_mafft_exe, gene_tree_seq, pwd_seq_file_1st_aln)
            for each_gene in SeqIO.parse(gene_tree_seq, 'fasta'):
                each_gene_genome = '_'.join(str(each_gene.id).split('_')[:-1])
                genome_subset.add(each_gene_genome)

        # run mafft
        os.system(cmd_mafft)

        # run fasttree
        # cmd_fasttree = '%s -quiet %s > %s 2>/dev/null' % (pwd_fasttree_exe, pwd_seq_file_2nd_aln, pwd_gene_tree_newick)
        cmd_fasttree = '%s -quiet %s > %s 2>/dev/null' % (pwd_fasttree_exe, pwd_seq_file_1st_aln, pwd_gene_tree_newick)
        os.system(cmd_fasttree)

        # Get species tree
        subset_tree(pwd_SCG_tree_all, genome_subset, pwd_species_tree_newick)

        # remove temp files
        os.remove(self_seq)
        os.remove(gene_tree_seq)
        if non_self_seq_num > 0:
            os.remove(non_self_seq)
            os.remove(blast_output)
            os.remove(blast_output_sorted)
            os.remove(gene_tree_seq_uniq)


def remove_hyphen_from_branch_length(tree_in, tree_out, tree_format):
    tree_in = Phylo.parse(tree_in, tree_format)
    Phylo.write(tree_in, tree_out, tree_format)


def Ranger_worker(argument_list):
    each_paired_tree = argument_list[0]
    pwd_ranger_inputs_folder = argument_list[1]
    pwd_tree_folder = argument_list[2]
    pwd_ranger_exe = argument_list[3]
    pwd_ranger_outputs_folder = argument_list[4]

    # define Ranger-DTL input file name
    each_paired_tree_concate = '___'.join(each_paired_tree)
    pwd_species_tree_newick                                 = '%s/%s_species_tree.newick'                               % (pwd_tree_folder, each_paired_tree_concate)
    pwd_gene_tree_newick                                    = '%s/%s_gene_tree.newick'                                  % (pwd_tree_folder, each_paired_tree_concate)
    pwd_species_tree_newick_no_hyphen_in_branch_length      = '%s/%s_species_tree_no_hyphen_in_branch_length.newick'    % (pwd_tree_folder, each_paired_tree_concate)
    pwd_gene_tree_newick_no_hyphen_in_branch_length         = '%s/%s_gene_tree_no_hyphen_in_branch_length.newick'       % (pwd_tree_folder, each_paired_tree_concate)

    # each_paired_tree_concate_short = '%s___%s' % (each_paired_tree[0].split('_')[-1], each_paired_tree[1].split('_')[-1])

    if (os.path.isfile(pwd_species_tree_newick) is True) and (os.path.isfile(pwd_gene_tree_newick) is True):

        ranger_inputs_file_name = each_paired_tree_concate + '.txt'
        ranger_outputs_file_name = each_paired_tree_concate + '_ranger_output.txt'

        pwd_ranger_inputs = '%s/%s' % (pwd_ranger_inputs_folder, ranger_inputs_file_name)
        pwd_ranger_outputs = '%s/%s' % (pwd_ranger_outputs_folder, ranger_outputs_file_name)


        # remove hyphen from branch length
        remove_hyphen_from_branch_length(pwd_species_tree_newick, pwd_species_tree_newick_no_hyphen_in_branch_length, 'newick')
        remove_hyphen_from_branch_length(pwd_gene_tree_newick, pwd_gene_tree_newick_no_hyphen_in_branch_length, 'newick')

        # read in species tree
        species_tree = Tree(pwd_species_tree_newick_no_hyphen_in_branch_length, format=0)
        species_tree.resolve_polytomy(recursive=True)  # solving multifurcations
        species_tree.convert_to_ultrametric()  # for dated mode

        # read in gene tree
        gene_tree = Tree(pwd_gene_tree_newick_no_hyphen_in_branch_length, format=0)
        gene_tree.resolve_polytomy(recursive=True)  # solving multifurcations

        ################################################################################################################

        # change species tree leaf name for Ranger-DTL2, replace "_" with "XXXXX", then, replace "." with "SSSSS", then replace "-" with "ZZZZZ"
        for each_st_leaf in species_tree:
            each_st_leaf_name = each_st_leaf.name

            # replace '_' with 'XXXXX'
            if '_' in each_st_leaf_name:
                each_st_leaf_name_no_Underline = 'XXAXX'.join(each_st_leaf_name.split('_'))
            else:
                each_st_leaf_name_no_Underline = each_st_leaf_name

            # replace '.' with 'SSSSS'
            if '.' in each_st_leaf_name_no_Underline:
                each_st_leaf_name_no_Underline_no_dot = 'SSASS'.join(each_st_leaf_name_no_Underline.split('.'))
            else:
                each_st_leaf_name_no_Underline_no_dot = each_st_leaf_name_no_Underline

            # replace '-' with 'ZZZZZ'
            if '-' in each_st_leaf_name_no_Underline_no_dot:
                each_st_leaf_name_no_Underline_no_dot_no_hyphen = 'ZZAZZ'.join(each_st_leaf_name_no_Underline_no_dot.split('-'))
            else:
                each_st_leaf_name_no_Underline_no_dot_no_hyphen = each_st_leaf_name_no_Underline_no_dot

            # rename species tree leaf name
            each_st_leaf.name = each_st_leaf_name_no_Underline_no_dot_no_hyphen

        # change gene tree leaf name for Ranger-DTL2, replace "_" with "XXXXX", then, replace "." with "SSSSS"
        for each_gt_leaf in gene_tree:
            each_gt_leaf_name = each_gt_leaf.name

            # replace '_' with 'XXXXX'
            if '_' in each_gt_leaf_name:
                each_gt_leaf_name_no_Underline = 'XXAXX'.join(each_gt_leaf_name.split('_')[:-1])
            else:
                each_gt_leaf_name_no_Underline = each_gt_leaf_name

            # replace '.' with 'SSSSS'
            if '.' in each_gt_leaf_name_no_Underline:
                each_gt_leaf_name_no_Underline_no_dot = 'SSASS'.join(each_gt_leaf_name_no_Underline.split('.'))
            else:
                each_gt_leaf_name_no_Underline_no_dot = each_gt_leaf_name_no_Underline

            # replace '-' with 'ZZZZZ'
            if '-' in each_gt_leaf_name_no_Underline_no_dot:
                each_gt_leaf_name_no_Underline_no_dot_no_hyphen = 'ZZAZZ'.join(each_gt_leaf_name_no_Underline_no_dot.split('-'))
            else:
                each_gt_leaf_name_no_Underline_no_dot_no_hyphen = each_gt_leaf_name_no_Underline_no_dot

            # rename gene tree leaf name
            each_gt_leaf.name = each_gt_leaf_name_no_Underline_no_dot_no_hyphen

        ################################################################################################################

        # write species tree and gene tree to Ranger-DTL input file
        ranger_inputs_file = open(pwd_ranger_inputs, 'w')

        # dated mode
        ranger_inputs_file.write('%s\n%s\n' % (species_tree.write(format=5), gene_tree.write(format=5)))
        ranger_inputs_file.close()

        # check if pwd_ranger_inputs is blank
        ranger_input_file_size = os.stat(pwd_ranger_inputs).st_size

        if ranger_input_file_size > 1:

            hyphen_detected = 'No'
            for line in open(pwd_ranger_inputs):
                if '-' in line:
                    hyphen_detected = 'Yes'

            # run Ranger-DTL
            if hyphen_detected == 'No':
                ranger_parameters = '-q -D 2 -T 3 -L 1'
                ranger_cmd = '%s %s -i %s -o %s' % (pwd_ranger_exe, ranger_parameters, pwd_ranger_inputs, pwd_ranger_outputs)
                os.system(ranger_cmd)


def extract_donor_recipient_sequences(pwd_combined_ffn, recipient_gene_list, donor_gene_list, pwd_recipient_gene_seq_ffn, pwd_recipient_gene_seq_faa, pwd_donor_gene_seq_ffn, pwd_donor_gene_seq_faa):

    pwd_recipient_gene_seq_ffn_handle = open(pwd_recipient_gene_seq_ffn, 'w')
    pwd_recipient_gene_seq_faa_handle = open(pwd_recipient_gene_seq_faa, 'w')
    pwd_donor_gene_seq_ffn_handle = open(pwd_donor_gene_seq_ffn, 'w')
    pwd_donor_gene_seq_faa_handle = open(pwd_donor_gene_seq_faa, 'w')

    for each_seq in SeqIO.parse(pwd_combined_ffn, 'fasta'):

        if str(each_seq.id) in recipient_gene_list:
            # write out nc sequences
            SeqIO.write(each_seq, pwd_recipient_gene_seq_ffn_handle, 'fasta')

            # write out aa sequences
            each_seq_aa = copy.deepcopy(each_seq)
            each_seq_aa.seq = each_seq_aa.seq.translate()
            SeqIO.write(each_seq_aa, pwd_recipient_gene_seq_faa_handle, 'fasta')

        if str(each_seq.id) in donor_gene_list:
            # write out nc sequences
            SeqIO.write(each_seq, pwd_donor_gene_seq_ffn_handle, 'fasta')

            # write out aa sequences
            each_seq_aa = copy.deepcopy(each_seq)
            each_seq_aa.seq = each_seq_aa.seq.translate()
            SeqIO.write(each_seq_aa, pwd_donor_gene_seq_faa_handle, 'fasta')

    pwd_recipient_gene_seq_ffn_handle.close()
    pwd_recipient_gene_seq_faa_handle.close()
    pwd_donor_gene_seq_ffn_handle.close()
    pwd_donor_gene_seq_faa_handle.close()


def Get_circlize_plot(multi_level_detection, output_prefix, pwd_candidates_file_PG_normal_txt, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, taxon_rank, taxon_rank_num, pwd_MetaCHIP_op_folder):

    rank_abbre_dict_plural      = {'d': 'domains', 'p': 'phyla', 'c': 'classes', 'o': 'orders', 'f': 'families', 'g': 'genera', 's': 'species', 'x': 'specified groups'}

    pwd_cir_plot_t1 =              '%s/%s_HGTs_among_%s_t1.txt'              % (pwd_MetaCHIP_op_folder, output_prefix, rank_abbre_dict_plural[taxon_rank])
    pwd_cir_plot_t1_sorted =       '%s/%s_HGTs_among_%s_t1_sorted.txt'       % (pwd_MetaCHIP_op_folder, output_prefix, rank_abbre_dict_plural[taxon_rank])
    pwd_cir_plot_t1_sorted_count = '%s/%s_HGTs_among_%s_sorted_count.txt'    % (pwd_MetaCHIP_op_folder, output_prefix, rank_abbre_dict_plural[taxon_rank])
    pwd_cir_plot_matrix_filename = '%s/%s_HGTs_among_%s.txt'                 % (pwd_MetaCHIP_op_folder, output_prefix, rank_abbre_dict_plural[taxon_rank])


    name2taxon_dict = {}
    transfers = []
    for each in open(pwd_candidates_file_PG_normal_txt):
        if not each.startswith('Gene_1'):
            each_split = each.strip().split('\t')
            Gene_1 = each_split[0]
            Gene_2 = each_split[1]
            Genome_1 = '_'.join(Gene_1.split('_')[:-1])
            Genome_2 = '_'.join(Gene_2.split('_')[:-1])

            if Genome_1 in genome_to_taxon_dict:
                Genome_1_taxon = '_'.join(genome_to_taxon_dict[Genome_1].split(' '))
            else:
                Genome_1_taxon = '%s_' % taxon_rank

            if Genome_2 in genome_to_taxon_dict:
                Genome_2_taxon = '_'.join(genome_to_taxon_dict[Genome_2].split(' '))
            else:
                Genome_2_taxon = '%s_' % taxon_rank

            Direction = each_split[5]
            if multi_level_detection == True:
                Direction = each_split[6]

            if '%)' in Direction:
                Direction = Direction.split('(')[0]

            if Genome_1 not in name2taxon_dict:
                name2taxon_dict[Genome_1] = Genome_1_taxon
            if Genome_2 not in name2taxon_dict:
                name2taxon_dict[Genome_2] = Genome_2_taxon
            transfers.append(Direction)


    tmp1 = open(pwd_cir_plot_t1, 'w')
    all_group_id = []
    for each_t in transfers:
        each_t_split = each_t.split('-->')
        donor = each_t_split[0]
        recipient = each_t_split[1]
        donor_id = name2taxon_dict[donor]
        recipient_id = name2taxon_dict[recipient]
        if donor_id not in all_group_id:
            all_group_id.append(donor_id)
        if recipient_id not in all_group_id:
            all_group_id.append(recipient_id)
        tmp1.write('%s,%s\n' % (donor_id, recipient_id))

    tmp1.close()

    os.system('cat %s | sort > %s' % (pwd_cir_plot_t1, pwd_cir_plot_t1_sorted))

    current_t = ''
    count = 0
    tmp2 = open(pwd_cir_plot_t1_sorted_count, 'w')
    for each_t2 in open(pwd_cir_plot_t1_sorted):
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
    for each_3 in open(pwd_cir_plot_t1_sorted_count):
        each_3_split = each_3.strip().split(',')
        key = '%s,%s' % (each_3_split[0], each_3_split[1])
        value = each_3_split[2]
        transfer_count[key] = value

    all_group_id = sorted(all_group_id)

    matrix_file = open(pwd_cir_plot_matrix_filename, 'w')
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
    if len(all_group_id) == 1:
        print('Too less group (1), plot skipped')
    elif 1 < len(all_group_id) <= 200:
        os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, pwd_cir_plot_matrix_filename, pwd_plot_circos))
    else:
        print('Too many groups (>200), plot skipped')

    # rm tmp files
    os.system('rm %s' % pwd_cir_plot_t1)
    os.system('rm %s' % pwd_cir_plot_t1_sorted)
    os.system('rm %s' % pwd_cir_plot_t1_sorted_count)


def Get_circlize_plot_customized_grouping(multi_level_detection, output_prefix, pwd_candidates_file_PG_normal_txt, genome_to_group_dict, circos_HGT_R, pwd_plot_circos, pwd_MetaCHIP_op_folder):

    pwd_cir_plot_t1 =              '%s/%s_cir_plot_t1.txt'              % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_t1_sorted =       '%s/%s_cir_plot_t1_sorted.txt'       % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_t1_sorted_count = '%s/%s_cir_plot_t1_sorted_count.txt' % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_matrix_filename = '%s/%s_cir_plot_matrix.csv'          % (pwd_MetaCHIP_op_folder, output_prefix)

    transfers = []
    for each in open(pwd_candidates_file_PG_normal_txt):
        if not each.startswith('Gene_1'):
            each_split = each.strip().split('\t')
            Gene_1 = each_split[0]
            Gene_2 = each_split[1]

            Direction = each_split[5]
            if multi_level_detection == True:
                Direction = each_split[6]

            if '%)' in Direction:
                Direction = Direction.split('(')[0]

            transfers.append(Direction)


    tmp1 = open(pwd_cir_plot_t1, 'w')
    all_group_id = []
    for each_t in transfers:
        each_t_split = each_t.split('-->')
        donor = each_t_split[0]
        recipient = each_t_split[1]
        donor_group = genome_to_group_dict[donor]
        recipient_group = genome_to_group_dict[recipient]
        if donor_group not in all_group_id:
            all_group_id.append(donor_group)
        if recipient_group not in all_group_id:
            all_group_id.append(recipient_group)
        tmp1.write('%s,%s\n' % (donor_group, recipient_group))

    tmp1.close()

    os.system('cat %s | sort > %s' % (pwd_cir_plot_t1, pwd_cir_plot_t1_sorted))

    current_t = ''
    count = 0
    tmp2 = open(pwd_cir_plot_t1_sorted_count, 'w')
    for each_t2 in open(pwd_cir_plot_t1_sorted):
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
    for each_3 in open(pwd_cir_plot_t1_sorted_count):
        each_3_split = each_3.strip().split(',')
        key = '%s,%s' % (each_3_split[0], each_3_split[1])
        value = each_3_split[2]
        transfer_count[key] = value

    all_group_id = sorted(all_group_id)

    matrix_file = open(pwd_cir_plot_matrix_filename, 'w')
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
    if len(all_group_id) > 1:
        os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, pwd_cir_plot_matrix_filename, pwd_plot_circos))

    # rm tmp files
    os.system('rm %s' % pwd_cir_plot_t1)
    os.system('rm %s' % pwd_cir_plot_t1_sorted)
    os.system('rm %s' % pwd_cir_plot_t1_sorted_count)


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def combine_PG_output(PG_output_file_list_with_path, output_prefix, detection_ranks, combined_PG_output_normal):

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


    combined_output_handle_normal = open(combined_PG_output_normal, 'w')

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
                    HGT_direction_freq = (concatenated_HGT_direction_list.count(HGT_direction)) * 100 / float(len(concatenated_HGT_direction_list))
                    if HGT_direction_freq > 50:
                        concatenated_HGT_direction = HGT_direction + '(' + str(float("{0:.2f}".format(HGT_direction_freq))) + '%)'

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

        if (HGT_end_match_dict[concatenated_HGT] == 'no') and (HGT_full_length_match_dict[concatenated_HGT] == 'no') and (concatenated_HGT_direction != 'NA') and (concatenated_HGT_direction != 'both'):
            combined_output_handle_normal.write(for_out)

    combined_output_handle_normal.close()


def BP(args, config_dict):

    output_folder =             args['o']
    output_prefix =             args['p']
    grouping_file =             args['g']
    grouping_levels =           args['r']
    cover_cutoff =              args['cov']
    align_len_cutoff =          args['al']
    flanking_length_kbp =       args['flk']
    identity_percentile =       args['ip']
    end_match_identity_cutoff = args['ei']
    num_threads =               args['t']
    No_Eb_Check =               args['NoEbCheck']
    keep_quiet =                args['quiet']
    keep_temp =                 args['tmp']


    # get path to current script
    flanking_length = flanking_length_kbp * 1000
    warnings.filterwarnings("ignore")
    rank_abbre_dict             = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species', 'x': 'specified group'}
    rank_abbre_dict_plural      = {'d': 'domains', 'p': 'phyla', 'c': 'classes', 'o': 'orders', 'f': 'families', 'g': 'genera', 's': 'species', 'x': 'specified groups'}
    rank_to_position_dict       = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}

    # check whether executables exist
    pwd_blastn_exe = config_dict['blastn']
    pwd_blastp_exe =    config_dict['blastp']
    pwd_mafft_exe =     config_dict['mafft']
    pwd_fasttree_exe =  config_dict['fasttree']
    circos_HGT_R =      config_dict['circos_HGT_R']
    pwd_ranger_exe = config_dict['ranger_linux']
    if platform.system() == 'Darwin':
        pwd_ranger_exe = config_dict['ranger_mac']

    program_list = [pwd_blastn_exe, pwd_ranger_exe, pwd_mafft_exe, pwd_fasttree_exe, pwd_blastp_exe]
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)
    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()

    #################################### find matched grouping file if not provided  ###################################

    if grouping_levels is None:
        grouping_levels = 'x'

    ################################### if provided ranks for PI and BP are different ##################################

    MetaCHIP_wd = '%s_MetaCHIP_wd' % output_prefix
    if output_folder is not None:
        MetaCHIP_wd = output_folder

    # get provided ranks for PI module
    pwd_prodigal_output_folder_re          = '%s/%s_*_prodigal_output' % (MetaCHIP_wd, output_prefix)
    pwd_prodigal_output_folder_detected    = [os.path.basename(file_name) for file_name in glob.glob(pwd_prodigal_output_folder_re)][0]
    specified_ranks_for_PI = pwd_prodigal_output_folder_detected.split('_prodigal_output')[0].split('_')[-1]


    specified_ranks_uniq_to_BP = []
    for each_BP_rank in grouping_levels:
        if each_BP_rank not in specified_ranks_for_PI:
            specified_ranks_uniq_to_BP.append(each_BP_rank)
    if len(specified_ranks_uniq_to_BP) >= 1:
        print('Not all provided ranks for BP module were found in PI ranks, program exited!')
        exit()

    if specified_ranks_for_PI != grouping_levels:
        prodigal_output_folder              = '%s_%s_prodigal_output'                       % (output_prefix, specified_ranks_for_PI)
        blast_db_folder                     = '%s_%s_blastdb'                               % (output_prefix, specified_ranks_for_PI)
        blast_cmd_file                      = '%s_%s_blastn_commands.txt'                   % (output_prefix, specified_ranks_for_PI)
        blast_result_folder                 = '%s_%s_blastn_results'                        % (output_prefix, specified_ranks_for_PI)
        combined_ffn_file                   = '%s_%s_combined_ffn.fasta'                    % (output_prefix, specified_ranks_for_PI)
        newick_tree_file                    = '%s_%s_SCG_tree.newick'                       % (output_prefix, specified_ranks_for_PI)
    else:
        prodigal_output_folder              = '%s_%s_prodigal_output'                       % (output_prefix, grouping_levels)
        blast_db_folder                     = '%s_%s_blastdb'                               % (output_prefix, grouping_levels)
        blast_cmd_file                      = '%s_%s_blastn_commands.txt'                   % (output_prefix, grouping_levels)
        blast_result_folder                 = '%s_%s_blastn_results'                        % (output_prefix, grouping_levels)
        combined_ffn_file                   = '%s_%s_combined_ffn.fasta'                    % (output_prefix, grouping_levels)
        newick_tree_file                    = '%s_%s_SCG_tree.newick'                       % (output_prefix, grouping_levels)

    ############################################# define file/folder names #############################################

    pwd_log_folder                      = '%s/%s_%s_log_files'                          % (MetaCHIP_wd, output_prefix, grouping_levels)
    pwd_flk_plots_folder                = '%s/%s_%s_flk_plots'                          % (MetaCHIP_wd, output_prefix, grouping_levels)
    pwd_log_file                        = '%s/%s_%s_BP_%s.log'                          % (pwd_log_folder, output_prefix, grouping_levels, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))
    pwd_ignored_taxonomic_rank_file     = '%s/ignored_taxonomic_rank.txt'               % MetaCHIP_wd
    blast_result_filtered_folder        = '%s_%s_blastn_results_filtered_al%sbp_cov%s'  % (output_prefix, grouping_levels, align_len_cutoff, cover_cutoff)
    pwd_prodigal_output_folder          = '%s/%s'                                       % (MetaCHIP_wd, prodigal_output_folder)
    pwd_blast_result_folder             = '%s/%s'                                       % (MetaCHIP_wd, blast_result_folder)
    pwd_blast_result_filtered_folder    = '%s/%s'                                       % (MetaCHIP_wd, blast_result_filtered_folder)
    pwd_combined_ffn_file               = '%s/%s/%s'                                    % (MetaCHIP_wd, blast_db_folder, combined_ffn_file)

    if os.path.isdir(pwd_log_folder) is False:
        os.mkdir(pwd_log_folder)

    ############################################### filter blastn results ##############################################

    # get ignored rank list
    ignored_rank_list = []
    if os.path.isfile(pwd_ignored_taxonomic_rank_file) is True:
        for ignored_rank in open(pwd_ignored_taxonomic_rank_file):
            ignored_rank_list.append(ignored_rank.strip())

    ignored_rank_num = 0
    for each_BP_rank in grouping_levels:
        if each_BP_rank in ignored_rank_list:
            ignored_rank_num += 1
    if ignored_rank_num == len(grouping_levels):
        print('All specified ranks were ignored, see report from the PI module.')
        exit()


    # check whether blast results exist
    blast_result_file_re                = '%s/*.tab' % pwd_blast_result_folder
    blast_result_file_list              = [os.path.basename(file_name) for file_name in glob.glob(blast_result_file_re)]
    blast_result_file_list_basename     = [i.split('_blastn.tab')[0] for i in blast_result_file_list]

    # exit if no blast result was found
    if len(blast_result_file_list) == 0:
        report_and_log(('No blastn reaults found in %s, please run all-vs-all blastn before the BP step!' % blast_result_folder), pwd_log_file, keep_quiet)
        exit()

    # exit if missing blast results were found
    ffn_file_re = '%s/*.ffn' % pwd_prodigal_output_folder
    ffn_file_list = [os.path.basename(file_name) for file_name in glob.glob(ffn_file_re)]
    ffn_file_list_basename = [i.split('.ffn')[0] for i in ffn_file_list]

    missing_blast_results = False
    for ffn_file in ffn_file_list_basename:
        if ffn_file not in blast_result_file_list_basename:
            report_and_log(('Blastn results for %s not found in folder %s!' % (ffn_file, blast_result_folder)), pwd_log_file, keep_quiet)
            missing_blast_results = True

    if missing_blast_results is True:
        report_and_log(('Program exited!'), pwd_log_file, keep_quiet)
        exit()

    # exit if the sizes of blast result files are all empty
    blast_result_file_list_empty = []
    blast_result_file_list_not_empty = []
    for blast_result_file in blast_result_file_list:
        pwd_blast_result_file = '%s/%s' % (pwd_blast_result_folder, blast_result_file)
        blast_result_file_size = os.stat(pwd_blast_result_file).st_size
        if blast_result_file_size == 0:
            blast_result_file_list_empty.append(blast_result_file)
        if blast_result_file_size > 0:
            blast_result_file_list_not_empty.append(blast_result_file)

    empty_file_percentage = float("{0:.2f}".format(len(blast_result_file_list_empty)*100/len(ffn_file_list)))

    if len(blast_result_file_list_empty) == len(ffn_file_list):
        report_and_log(('All files in %s are empty, program exited' % blast_result_folder), pwd_log_file, keep_quiet)
        exit()
    elif empty_file_percentage > 20:
        report_and_log(('%s%s of the blastn results files in %s are empty, program exited!' % (empty_file_percentage, '%', blast_result_folder)), pwd_log_file, keep_quiet)
        exit()
    elif 0 < empty_file_percentage <= 20:
        for empty_file in blast_result_file_list_empty:
            report_and_log(('%s is empty, HGT detection will be performed based on the rest files.' % empty_file), pwd_log_file, keep_quiet)

    # check whether filtered blast results already exist
    filtered_blast_results_found = False
    if os.path.isdir(pwd_blast_result_filtered_folder) is True:
        blast_result_filtered_file_re = '%s/*_filtered.tab' % pwd_blast_result_filtered_folder
        blast_result_filtered_file_list = glob.glob(blast_result_filtered_file_re)

        if len(blast_result_filtered_file_list) == len(ffn_file_list):
            filtered_blast_results_found = True

    if filtered_blast_results_found is True:
        report_and_log(('Filtered blastn results detected, filtration step skipped.'), pwd_log_file, keep_quiet)
    else:
        report_and_log(('Filtering blastn results with the following criteria: Query genome != Subject genome, Alignment length >= %sbp and coverage >= %s%s.' % (align_len_cutoff, cover_cutoff, '%')), pwd_log_file, keep_quiet)

        # create folder
        force_create_folder(pwd_blast_result_filtered_folder)

        # filter blastn results with multiprocessing
        list_for_multiple_arguments_filter_blast_results = []
        for blast_result_file in blast_result_file_list:
            pwd_blast_result_file = '%s/%s' % (pwd_blast_result_folder, blast_result_file)
            blast_result_filtered_file = '%s_filtered.tab' % ('.'.join(blast_result_file.split('.')[:-1]))
            pwd_blast_result_filtered_file = '%s/%s' % (pwd_blast_result_filtered_folder, blast_result_filtered_file)
            list_for_multiple_arguments_filter_blast_results.append([pwd_blast_result_file, align_len_cutoff, cover_cutoff, pwd_blast_result_filtered_file])

        # filter_blast_results with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(filter_blast_results_worker, list_for_multiple_arguments_filter_blast_results)
        pool.close()
        pool.join()

    ############################################## perform HGT detection ###############################################

    # perform BM and PG approaches at all specified taxonomic ranks
    ploted_flk_plot_list = []
    customized_group_num = 0
    for grouping_level in grouping_levels:

        if grouping_level not in ignored_rank_list:

            # read in grouping information
            pwd_grouping_file = ''
            group_num = 0
            if grouping_file is None:
                grouping_file_re = '%s/%s_grouping_%s*.txt' % (MetaCHIP_wd, output_prefix, grouping_level)
                grouping_file_list = [os.path.basename(file_name) for file_name in glob.glob(grouping_file_re)]
                if len(grouping_file_list) == 1:
                    detected_grouping_file = grouping_file_list[0]
                    pwd_grouping_file = '%s/%s' % (MetaCHIP_wd, detected_grouping_file)
                    group_num = get_number_of_group(pwd_grouping_file)
                    report_and_log(('Detect HGT among %s: input genomes were clustered into %s %s.' % (rank_abbre_dict_plural[grouping_level], group_num, rank_abbre_dict_plural[grouping_level])), pwd_log_file, keep_quiet)

                elif len(grouping_file_list) == 0:
                    report_and_log(('Detect HGT among %s: no grouping file at %s level found, program exited.' % (rank_abbre_dict_plural[grouping_level], rank_abbre_dict[grouping_level])), pwd_log_file, keep_quiet)
                    exit()

                else:
                    report_and_log(('Detect HGT among %s: multiple grouping file at %s level found, program exited.' % (rank_abbre_dict_plural[grouping_level], rank_abbre_dict[grouping_level])), pwd_log_file, keep_quiet)
                    exit()

            else:  # with provided grouping file
                pwd_grouping_file = grouping_file
                group_num = get_number_of_group(pwd_grouping_file)

                customized_group_num = group_num


            ############################################# define file/folder names #############################################

            MetaCHIP_op_folder                          = '%s_HGT_ip%s_al%sbp_c%s_ei%s_f%skbp_%s%s'         % (output_prefix, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff),str(end_match_identity_cutoff), flanking_length_kbp, grouping_level, group_num)
            blast_result_filtered_folder_g2g            = '%s_%s%s_1_blastn_results_filtered_g2g'           % (output_prefix, grouping_level, group_num)
            grouping_file_with_id_filename              = '%s_%s%s_grouping_with_id.txt'                    % (output_prefix, grouping_level, group_num)
            qual_idens_file_gg                          = '%s_%s%s_qualified_iden_gg.txt'                   % (output_prefix, grouping_level, group_num)
            qual_idens_file_gg_sorted                   = '%s_%s%s_qualified_iden_gg_sorted.txt'            % (output_prefix, grouping_level, group_num)
            blast_result_filtered_folder_with_group     = '%s_%s%s_2_blastn_results_filtered_with_group'    % (output_prefix, grouping_level, group_num)
            blast_result_filtered_folder_in_one_line    = '%s_%s%s_3_blastn_results_filtered_in_one_line'   % (output_prefix, grouping_level, group_num)
            op_candidates_only_gene_folder_name         = '%s_%s%s_4_HGTs_only_id'                          % (output_prefix, grouping_level, group_num)
            op_candidates_with_group_folder_name        = '%s_%s%s_5_HGTs_with_group'                       % (output_prefix, grouping_level, group_num)
            op_candidates_with_group_file_name          = '%s_%s%s_HGTs_with_group.txt'                     % (output_prefix, grouping_level, group_num)
            op_candidates_only_gene_file_name           = '%s_%s%s_HGTs_only_id.txt'                        % (output_prefix, grouping_level, group_num)
            op_candidates_only_gene_file_name_uniq      = '%s_%s%s_HGTs_uniq.txt'                           % (output_prefix, grouping_level, group_num)
            op_act_folder_name                          = '%s_%s%s_Flanking_region_plots'                   % (output_prefix, grouping_level, group_num)
            normal_folder_name                          = '1_Plots_normal'
            end_match_folder_name                       = '2_Plots_end_match'
            full_length_match_folder_name               = '3_Plots_full_length_match'
            op_candidates_BM                            = '%s_%s%s_HGTs_BM.txt'                             % (output_prefix, grouping_level, group_num)
            subjects_in_one_line_filename               = '%s_%s%s_subjects_in_one_line.txt'                % (output_prefix, grouping_level, group_num)
            HGT_query_to_subjects_filename              = '%s_%s%s_HGT_query_to_subjects.txt'               % (output_prefix, grouping_level, group_num)
            op_candidates_seq_nc                        = '%s_%s%s_HGTs_BM_nc.fasta'                        % (output_prefix, grouping_level, group_num)
            combined_faa_file                           = '%s_%s_combined_faa.fasta'                        % (output_prefix, grouping_levels)
            combined_faa_file_subset                    = '%s_%s%s_combined_subset.faa'                     % (output_prefix, grouping_level, group_num)
            tree_folder                                 = '%s_%s%s_PG_tree_folder'                          % (output_prefix, grouping_level, group_num)
            ranger_inputs_folder_name                   = '%s_%s%s_PG_Ranger_input'                         % (output_prefix, grouping_level, group_num)
            ranger_outputs_folder_name                  = '%s_%s%s_PG_Ranger_output'                        % (output_prefix, grouping_level, group_num)
            candidates_file_name                        = '%s_%s%s_HGTs_BM.txt'                             % (output_prefix, grouping_level, group_num)
            candidates_seq_file_name                    = '%s_%s%s_HGTs_BM_nc.fasta'                        % (output_prefix, grouping_level, group_num)
            candidates_file_name_ET                     = '%s_%s%s_HGTs_PG.txt'                             % (output_prefix, grouping_level, group_num)
            flanking_region_plot_folder_name            = '%s_%s%s_Flanking_region_plots'                   % (output_prefix, grouping_level, group_num)
            grouping_file_with_id_filename              = '%s_%s%s_grouping_with_id.txt'                    % (output_prefix, grouping_level, group_num)
            HGT_query_to_subjects_filename              = '%s_%s%s_HGT_query_to_subjects.txt'               % (output_prefix, grouping_level, group_num)

            pwd_MetaCHIP_op_folder                          = '%s/%s'       % (MetaCHIP_wd, MetaCHIP_op_folder)
            pwd_blast_result_filtered_folder_g2g            = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, blast_result_filtered_folder_g2g)
            pwd_grouping_file_with_id                       = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, grouping_file_with_id_filename)
            pwd_qual_iden_file_gg                           = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, qual_idens_file_gg)
            pwd_qual_iden_file_gg_sorted                    = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, qual_idens_file_gg_sorted)
            pwd_blast_result_filtered_folder_with_group     = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, blast_result_filtered_folder_with_group)
            pwd_blast_result_filtered_folder_in_one_line    = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, blast_result_filtered_folder_in_one_line)
            pwd_op_candidates_with_group_folder             = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_with_group_folder_name)
            pwd_op_candidates_only_gene_folder              = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_only_gene_folder_name)
            pwd_op_candidates_with_group_file               = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_with_group_file_name)
            pwd_op_candidates_only_gene_file                = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_only_gene_file_name)
            pwd_op_candidates_only_gene_file_uniq           = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_only_gene_file_name_uniq)
            pwd_op_act_folder                               = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name)
            pwd_normal_folder                               = '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name, normal_folder_name)
            pwd_end_match_folder                            = '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name, end_match_folder_name)
            pwd_full_length_match_folder                    = '%s/%s/%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder, op_act_folder_name, full_length_match_folder_name)
            pwd_op_candidates_BM                            = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_BM)
            pwd_subjects_in_one_line                        = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, subjects_in_one_line_filename)
            pwd_HGT_query_to_subjects_file                  = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, HGT_query_to_subjects_filename)
            pwd_op_candidates_seq_nc                        = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, op_candidates_seq_nc)
            pwd_candidates_file                             = '%s/%s'       % (pwd_MetaCHIP_op_folder, candidates_file_name)
            pwd_candidates_seq_file                         = '%s/%s'       % (pwd_MetaCHIP_op_folder, candidates_seq_file_name)
            pwd_candidates_file_ET                          = '%s/%s'       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET)
            pwd_combined_faa_file                           = '%s/%s'       % (pwd_MetaCHIP_op_folder, combined_faa_file)
            pwd_ranger_inputs_folder                        = '%s/%s'       % (pwd_MetaCHIP_op_folder, ranger_inputs_folder_name)
            pwd_ranger_outputs_folder                       = '%s/%s'       % (pwd_MetaCHIP_op_folder, ranger_outputs_folder_name)
            pwd_tree_folder                                 = '%s/%s'       % (pwd_MetaCHIP_op_folder, tree_folder)
            pwd_combined_faa_file_subset                    = '%s/%s'       % (pwd_MetaCHIP_op_folder, combined_faa_file_subset)
            pwd_grouping_file_with_id                       = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, grouping_file_with_id_filename)
            pwd_HGT_query_to_subjects_file                  = '%s/%s/%s'    % (MetaCHIP_wd, MetaCHIP_op_folder, HGT_query_to_subjects_filename)
            pwd_newick_tree_file                            = '%s/%s'       % (MetaCHIP_wd, newick_tree_file)

            force_create_folder(pwd_MetaCHIP_op_folder)


            ####################################################################################################################

            # index grouping file
            index_grouping_file(pwd_grouping_file, pwd_grouping_file_with_id)

            # create genome_group_dict and genome_list
            qualified_genome_list = []
            name_to_group_number_dict = {}
            name_to_group_dict = {}
            for each_bin in open(pwd_grouping_file_with_id):
                each_bin_split = each_bin.strip().split(',')
                bin_name = each_bin_split[1]
                bin_group_number = each_bin_split[0]
                bin_group = bin_group_number.split('_')[0]
                qualified_genome_list.append(bin_name)
                name_to_group_number_dict[bin_name] = bin_group_number
                name_to_group_dict[bin_name] = bin_group


            ############################################################################################################
            ################################## perform HGT detection with BM approach ##################################
            ############################################################################################################

            report_and_log(('Detect HGT among %s: Best-match approach.' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)

            ########################################### get group-to-group identities ##########################################

            os.mkdir(pwd_blast_result_filtered_folder_g2g)

            report_and_log(('Detect HGT among %s: get group-to-group identities with %s cores.' % (rank_abbre_dict_plural[grouping_level], num_threads)), pwd_log_file, keep_quiet)

            # create folder
            force_create_folder(pwd_blast_result_filtered_folder_g2g)

            # get file list
            blast_result_filtered_file_re = '%s/*_filtered.tab' % pwd_blast_result_filtered_folder
            blast_result_filtered_file_list = [os.path.basename(file_name) for file_name in glob.glob(blast_result_filtered_file_re)]

            list_for_multiple_arguments_get_g2g_identities = []
            for filtered_blast_result in blast_result_filtered_file_list:
                genome_id = filtered_blast_result.split('_blastn_filtered.tab')[0]
                if genome_id in qualified_genome_list:
                    pwd_filtered_blast_result = '%s/%s' % (pwd_blast_result_filtered_folder, filtered_blast_result)
                    pwd_filtered_blast_result_g2g = '%s/%s_g2g.tab' % (pwd_blast_result_filtered_folder_g2g, filtered_blast_result.split('.tab')[0])
                    list_for_multiple_arguments_get_g2g_identities.append([pwd_filtered_blast_result, qualified_genome_list, name_to_group_number_dict, pwd_filtered_blast_result_g2g])

            # get group-to-group identities with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(get_g2g_identities_worker, list_for_multiple_arguments_get_g2g_identities)
            pool.close()
            pool.join()

            # combine g2g_files and sort it
            os.system('cat %s/*_g2g.tab > %s' % (pwd_blast_result_filtered_folder_g2g, pwd_qual_iden_file_gg))
            os.system('cat %s | sort > %s' % (pwd_qual_iden_file_gg, pwd_qual_iden_file_gg_sorted))


            ###########################################  get cutoff according to specified percentile ###########################################

            # get identities for each group pair and get group_pair to identity dict
            current_group_pair_name = ''
            current_group_pair_identities = []
            group_pair_iden_cutoff_dict = {}
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
                        current_group_pair_identities_array = np.array(current_group_pair_identities)
                        current_group_pair_identity_cut_off = np.percentile(current_group_pair_identities_array,identity_percentile)
                        current_group_pair_identity_cut_off = float("{0:.2f}".format(current_group_pair_identity_cut_off))
                        current_group_pair_name_split = current_group_pair_name.split('_')
                        current_group_pair_name_swapped = '%s_%s' % (current_group_pair_name_split[1], current_group_pair_name_split[0])
                        group_pair_iden_cutoff_dict[current_group_pair_name] = current_group_pair_identity_cut_off
                        group_pair_iden_cutoff_dict[current_group_pair_name_swapped] = current_group_pair_identity_cut_off

                        # reset
                        current_group_pair_name = group_pair
                        current_group_pair_identities = []
                        current_group_pair_identities.append(identity)

            # for the last group
            current_group_pair_identities_array = np.array(current_group_pair_identities)
            current_group_pair_identity_cut_off = np.percentile(current_group_pair_identities_array, identity_percentile)
            current_group_pair_identity_cut_off = float("{0:.2f}".format(current_group_pair_identity_cut_off))
            current_group_pair_name_split = current_group_pair_name.split('_')
            current_group_pair_name_swapped = '%s_%s' % (current_group_pair_name_split[1], current_group_pair_name_split[0])
            group_pair_iden_cutoff_dict[current_group_pair_name] = current_group_pair_identity_cut_off
            group_pair_iden_cutoff_dict[current_group_pair_name_swapped] = current_group_pair_identity_cut_off


            ############################### add group to blast hits and put subjects in one line ###############################

            report_and_log(('Detect HGT among %s: analyzing Blast hits with %s cores.' % (rank_abbre_dict_plural[grouping_level], num_threads)), pwd_log_file, keep_quiet)

            # create folder
            force_create_folder(pwd_blast_result_filtered_folder_with_group)
            force_create_folder(pwd_blast_result_filtered_folder_in_one_line)
            force_create_folder(pwd_op_candidates_with_group_folder)
            force_create_folder(pwd_op_candidates_only_gene_folder)

            list_for_multiple_arguments_get_HGT = []
            for filtered_blast_result in blast_result_filtered_file_list:
                genome_id = filtered_blast_result.split('_blastn_filtered.tab')[0]
                if genome_id in qualified_genome_list:
                    pwd_filtered_blast_result               = '%s/%s'                           % (pwd_blast_result_filtered_folder, filtered_blast_result)
                    pwd_filtered_blast_result_with_group    = '%s/%s_with_group_sorted.tab'     % (pwd_blast_result_filtered_folder_with_group, filtered_blast_result.split('.tab')[0])
                    pwd_filtered_blast_result_in_one_line   = '%s/%s_subjects_in_one_line.tab'  % (pwd_blast_result_filtered_folder_in_one_line, filtered_blast_result.split('.tab')[0])
                    pwd_hgt_candidates_with_group           = '%s/%s_HGTs_with_group.txt'       % (pwd_op_candidates_with_group_folder, genome_id)
                    pwd_hgt_candidates_only_gene            = '%s/%s_HGTs_only_gene.txt'        % (pwd_op_candidates_only_gene_folder, genome_id)

                    list_for_multiple_arguments_get_HGT.append([pwd_filtered_blast_result,
                                                                name_to_group_number_dict,
                                                                pwd_filtered_blast_result_with_group,
                                                                pwd_filtered_blast_result_in_one_line,
                                                                pwd_hgt_candidates_with_group,
                                                                pwd_hgt_candidates_only_gene,
                                                                group_pair_iden_cutoff_dict,
                                                                qualified_genome_list])

            # add group to blast hits with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(get_HGT_worker, list_for_multiple_arguments_get_HGT)
            pool.close()
            pool.join()


            ################################ remove bidirection and add identity to output file ################################

            # combine pwd_op_candidates_with_group_file and pwd_op_candidates_only_gene_file
            os.system('cat %s/*_HGTs_with_group.txt > %s' % (pwd_op_candidates_with_group_folder, pwd_op_candidates_with_group_file))
            os.system('cat %s/*_HGTs_only_gene.txt > %s' % (pwd_op_candidates_only_gene_folder, pwd_op_candidates_only_gene_file))

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

            report_and_log(('Detect HGT among %s: plotting flanking regions with %s cores.' % (rank_abbre_dict_plural[grouping_level], num_threads)), pwd_log_file, keep_quiet)

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
                list_for_multiple_arguments_flanking_regions.append([match,
                                                                     pwd_prodigal_output_folder,
                                                                     flanking_length,
                                                                     align_len_cutoff,
                                                                     name_to_group_number_dict,
                                                                     pwd_op_act_folder,
                                                                     pwd_normal_folder,
                                                                     pwd_end_match_folder,
                                                                     pwd_full_length_match_folder,
                                                                     pwd_blastn_exe,
                                                                     keep_temp,
                                                                     candidates_2_contig_match_category_dict_mp,
                                                                     end_match_identity_cutoff,
                                                                     No_Eb_Check])

            pool_flanking_regions = mp.Pool(processes=num_threads)
            pool_flanking_regions.map(get_gbk_blast_act2, list_for_multiple_arguments_flanking_regions)
            pool_flanking_regions.close()
            pool_flanking_regions.join()

            # remove temporary folder
            if keep_temp == 0:
                act_file_re = '%s/*___*/*' % pwd_op_act_folder
                act_folder_re = '%s/*___*' % pwd_op_act_folder
                rm_folder_file(act_file_re)
                rm_folder_file(act_folder_re)

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
                BM_output_file_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                recipient_gene, donor_gene, recipient_genome_group, donor_genome_group, identity, end_match_value,
                full_length_match_value))
            BM_output_file_handle.close()


            ####################################### export gene clusters for PG approach #######################################

            os.system('cat %s/*_in_one_line.tab > %s' % (pwd_blast_result_filtered_folder_in_one_line, pwd_subjects_in_one_line))
            export_HGT_query_to_subjects(pwd_op_candidates_BM, pwd_subjects_in_one_line, pwd_HGT_query_to_subjects_file)


            ################################### export nc and aa sequence of predicted HGTs ####################################

            report_and_log(('Detect HGT among %s: get sequences of BM predicted HGT candidates.' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)

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

            # report
            report_and_log(('Detect HGT among %s: done for BM approach!' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)


            ############################################################################################################
            ################################## perform HGT detection with PG approach ##################################
            ############################################################################################################

            ###################################### store ortholog information into dictionary ######################################

            # create folders
            force_create_folder(pwd_tree_folder)

            # get list of match pair list
            candidates_list = []
            candidates_list_genes = set()
            for match_group in open(pwd_candidates_file):
                if not match_group.startswith('Gene_1'):
                    match_group_split = match_group.strip().split('\t')
                    end_match = match_group_split[5]
                    full_length_match = match_group_split[6]
                    if (end_match == 'no') and (full_length_match == 'no'):
                        candidates_list.append(match_group_split[:2])
                        candidates_list_genes.add(match_group_split[0])
                        candidates_list_genes.add(match_group_split[1])

            if candidates_list == []:
                report_and_log(('Detect HGT among %s: no HGT detected by BM approach!' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet=False)

            # for report and log
            report_and_log(('Detect HGT among %s: get gene/genome members in gene/species tree for BM predicted HGT candidates.' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)

            # get bin_record_list and genome name list
            bin_record_list = []
            genome_name_list = []
            name_to_group_dict = {}
            name_to_group_number_dict = {}
            name_to_group_number_without_underscore_dict = {}
            bin_group_list = []
            bin_group_without_underscore_list = []
            for each_bin in open(pwd_grouping_file_with_id):
                each_bin_split = each_bin.strip().split(',')
                bin_group = each_bin_split[0]
                bin_group_without_underscore = bin_group.split('_')[0] + bin_group.split('_')[1]
                bin_name = each_bin_split[1]
                name_to_group_dict[bin_name] = bin_group.split('_')[0]
                name_to_group_number_dict[bin_name] = bin_group
                name_to_group_number_without_underscore_dict[bin_name] = bin_group_without_underscore
                bin_record = BinRecord(bin_name, bin_group, bin_group_without_underscore)
                bin_record_list.append(bin_record)
                genome_name_list.append(bin_name)
                bin_group_list.append(bin_group)
                bin_group_without_underscore_list.append(bin_group_without_underscore)


            ###################################################### Get dicts #######################################################

            # get HGT_query_to_subjects dict
            HGT_query_to_subjects_dict = {}
            gene_id_overall = set()
            for each_candidate in open(pwd_HGT_query_to_subjects_file):
                each_candidate_split = each_candidate.strip().split('\t')
                query = each_candidate_split[0]
                subjects = each_candidate_split[1].split(',')
                if query in candidates_list_genes:
                    HGT_query_to_subjects_dict[query] = subjects
                    for each_subject in subjects:
                        gene_id_overall.add(each_subject)


            ################################# Prepare subset of faa_file for building gene tree ####################################

            # for report and log
            report_and_log(('Detect HGT among %s: prepare %s subset.' % (rank_abbre_dict_plural[grouping_level], combined_faa_file)), pwd_log_file, keep_quiet)

            # uniq gene id list
            gene_id_uniq_set = set()
            for each_gene in gene_id_overall:
                gene_id_uniq_set.add(each_gene)

            # get combined_faa_file file
            os.system('cat %s/*.faa > %s' % (pwd_prodigal_output_folder, pwd_combined_faa_file))

            # prepare combined_ffn file subset to speed up
            pwd_combined_faa_file_subset_handle = open(pwd_combined_faa_file_subset, 'w')
            for each_gene in SeqIO.parse(pwd_combined_faa_file, 'fasta'):
                if each_gene.id in gene_id_uniq_set:
                    pwd_combined_faa_file_subset_handle.write('>%s\n' % each_gene.id)
                    pwd_combined_faa_file_subset_handle.write('%s\n' % str(each_gene.seq))
            pwd_combined_faa_file_subset_handle.close()

            ################################## get gene tree and SCG tree subset ##################################

            # for report and log
            report_and_log(('Detect HGT among %s: get species and gene tree for %s BM approach identified HGTs.' % (rank_abbre_dict_plural[grouping_level], len(candidates_list))), pwd_log_file, keep_quiet)

            # put multiple arguments in list
            list_for_multiple_arguments_extract_gene_tree_seq = []
            for each_to_extract in candidates_list:
                list_for_multiple_arguments_extract_gene_tree_seq.append([each_to_extract,
                                                                          pwd_tree_folder,
                                                                          pwd_combined_faa_file_subset,
                                                                          pwd_blastp_exe,
                                                                          pwd_mafft_exe,
                                                                          pwd_fasttree_exe,
                                                                          name_to_group_dict,
                                                                          genome_name_list,
                                                                          HGT_query_to_subjects_dict,
                                                                          pwd_newick_tree_file])
            pool = mp.Pool(processes=num_threads)
            pool.map(extract_gene_tree_seq_worker, list_for_multiple_arguments_extract_gene_tree_seq)
            pool.close()
            pool.join()

            ##################################################### Run Ranger-DTL ###################################################

            # prepare folders
            force_create_folder(pwd_ranger_inputs_folder)
            force_create_folder(pwd_ranger_outputs_folder)

            # for report and log
            report_and_log(('Detect HGT among %s: running Ranger-DTL2.'% rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)

            # put multiple arguments in list
            list_for_multiple_arguments_Ranger = []
            for each_paired_tree in candidates_list:
                list_for_multiple_arguments_Ranger.append(
                    [each_paired_tree, pwd_ranger_inputs_folder, pwd_tree_folder, pwd_ranger_exe,
                     pwd_ranger_outputs_folder])

            pool = mp.Pool(processes=num_threads)
            pool.map(Ranger_worker, list_for_multiple_arguments_Ranger)
            pool.close()
            pool.join()

            ########################################### parse Ranger-DTL prediction result #########################################

            # for report and log
            report_and_log(('Detect HGT among %s: parsing Ranger-DTL2 outputs.' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)

            candidate_2_predictions_dict = {}
            candidate_2_possible_direction_dict = {}
            for each_ranger_prediction in candidates_list:
                each_ranger_prediction_concate = '___'.join(each_ranger_prediction)
                ranger_out_file_name = each_ranger_prediction_concate + '_ranger_output.txt'
                pwd_ranger_result = '%s/%s' % (pwd_ranger_outputs_folder, ranger_out_file_name)

                if os.path.isfile(pwd_ranger_result) == True:

                    # parse prediction result
                    predicted_transfers = []
                    for each_line in open(pwd_ranger_result):
                        if 'Transfer' in each_line:
                            if not each_line.startswith('The minimum reconciliation cost'):
                                mapping = each_line.strip().split(':')[1].split(',')[1]
                                recipient = each_line.strip().split(':')[1].split(',')[2]
                                donor_p = mapping.split('-->')[1][1:]
                                donor_p = '_'.join(donor_p.split('XXAXX'))
                                donor_p = '.'.join(donor_p.split('SSASS'))
                                donor_p = '-'.join(donor_p.split('ZZAZZ'))

                                recipient_p = recipient.split('-->')[1][1:]
                                recipient_p = '_'.join(recipient_p.split('XXAXX'))
                                recipient_p = '.'.join(recipient_p.split('SSASS'))
                                recipient_p = '-'.join(recipient_p.split('ZZAZZ'))
                                predicted_transfer = donor_p + '-->' + recipient_p
                                predicted_transfers.append(predicted_transfer)

                    candidate_2_predictions_dict[each_ranger_prediction_concate] = predicted_transfers

                    # get two possible transfer situation
                    candidate_split_gene = each_ranger_prediction_concate.split('___')
                    candidate_split_gene_only_genome = []
                    for each_candidate in candidate_split_gene:
                        each_candidate_genome = '_'.join(each_candidate.split('_')[:-1])
                        candidate_split_gene_only_genome.append(each_candidate_genome)

                    possible_hgt_1 = '%s-->%s' % (
                    candidate_split_gene_only_genome[0], candidate_split_gene_only_genome[1])
                    possible_hgt_2 = '%s-->%s' % (
                    candidate_split_gene_only_genome[1], candidate_split_gene_only_genome[0])
                    possible_hgts = [possible_hgt_1, possible_hgt_2]
                    candidate_2_possible_direction_dict[each_ranger_prediction_concate] = possible_hgts


            #################################################### combine results ###################################################

            # for report and log
            report_and_log(('Detect HGT among %s: add Ranger-DTL predicted direction to %s' % (rank_abbre_dict_plural[grouping_level], candidates_file_name_ET)), pwd_log_file, keep_quiet)

            # add results to output file of best blast match approach
            combined_output_handle = open(pwd_candidates_file_ET, 'w')
            # combined_output_validated_handle = open(pwd_candidates_file_ET_validated, 'w')
            combined_output_validated_header = 'Gene_1\tGene_2\tGene_1_group\tGene_2_group\tIdentity\tend_match\tfull_length_match\tDirection\n'
            # combined_output_validated_handle.write(combined_output_validated_header)
            combined_output_handle.write(combined_output_validated_header)
            validated_candidate_list = []
            for match_group in open(pwd_candidates_file):
                if not match_group.startswith('Gene_1'):
                    match_group_split = match_group.strip().split('\t')
                    recipient_gene = match_group_split[0]
                    donor_gene = match_group_split[1]
                    recipient_genome_id = match_group_split[2]
                    donor_genome_id = match_group_split[3]
                    identity = match_group_split[4]
                    end_break = match_group_split[5]
                    Ctg_align = match_group_split[6]
                    concatenated = '%s___%s' % (recipient_gene, donor_gene)
                    possible_direction = []
                    if concatenated in candidate_2_possible_direction_dict:
                        possible_direction = candidate_2_possible_direction_dict[concatenated]

                    validated_prediction = 'NA'
                    if concatenated in candidate_2_predictions_dict:
                        for each_prediction in candidate_2_predictions_dict[concatenated]:
                            if each_prediction in possible_direction:
                                validated_prediction = each_prediction

                    if (Ctg_align == 'no') and (end_break == 'no') and (validated_prediction != 'NA'):
                        if recipient_gene not in validated_candidate_list:
                            validated_candidate_list.append(recipient_gene)
                        if donor_gene not in validated_candidate_list:
                            validated_candidate_list.append(donor_gene)
                        # combined_output_validated_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_id, donor_genome_id, identity, end_break, Ctg_align, validated_prediction))
                    combined_output_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_id, donor_genome_id, identity, end_break, Ctg_align, validated_prediction))
            combined_output_handle.close()


            ################################################### remove tmp files ###################################################

            #report_and_log(('Detect HGT among %s: done for Phylogenetic approach!' % rank_abbre_dict_plural[grouping_level]), pwd_log_file, keep_quiet)
            if keep_temp is False:
                shutil.rmtree(pwd_blast_result_filtered_folder_g2g)
                shutil.rmtree(pwd_blast_result_filtered_folder_with_group)
                shutil.rmtree(pwd_blast_result_filtered_folder_in_one_line)
                shutil.rmtree(pwd_op_candidates_only_gene_folder)
                shutil.rmtree(pwd_op_candidates_with_group_folder)
                shutil.rmtree(pwd_ranger_inputs_folder)
                shutil.rmtree(pwd_ranger_outputs_folder)
                shutil.rmtree(pwd_tree_folder)

                os.remove(pwd_HGT_query_to_subjects_file)
                os.remove(pwd_subjects_in_one_line)
                os.remove(pwd_qual_iden_file_gg)
                os.remove(pwd_qual_iden_file_gg_sorted)
                os.remove(pwd_op_candidates_with_group_file)
                os.remove(pwd_op_candidates_only_gene_file)
                os.remove(pwd_op_candidates_only_gene_file_uniq)
                os.remove(pwd_grouping_file_with_id)
                os.remove(pwd_combined_faa_file_subset)


    ####################################################################################################################
    ########################################### combine BM and PG predictions ##########################################
    ####################################################################################################################

    if grouping_file is not None:

        multi_level_detection = False
        pwd_MetaCHIP_op_folder_re = '%s/%s_HGT_ip%s_al%sbp_c%s_ei%s_f%skbp_x*' % (MetaCHIP_wd, output_prefix, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)

        MetaCHIP_op_folder = ''
        if len([os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)]) == 1:
            MetaCHIP_op_folder = [os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)][0]
        else:
            for op_folder in [os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)]:
                group_number = int(op_folder.split('_')[-1][1:])
                if group_number == customized_group_num:
                    MetaCHIP_op_folder = op_folder

        group_num = int(MetaCHIP_op_folder.split('_')[-1][1:])

        pwd_MetaCHIP_op_folder =        '%s/%s' % (MetaCHIP_wd, MetaCHIP_op_folder)
        pwd_detected_HGT_PG_txt =       '%s/%s_x%s_HGTs_PG.txt'                      % (pwd_MetaCHIP_op_folder, output_prefix, group_num)
        pwd_flanking_plot_folder =      '%s/%s_x%s_Flanking_region_plots'            % (pwd_MetaCHIP_op_folder, output_prefix, group_num)
        pwd_detected_HGT_txt =          '%s/%s_x_detected_HGTs.txt'                   % (pwd_MetaCHIP_op_folder, output_prefix)
        pwd_recipient_gene_seq_ffn =    '%s/%s_x_detected_HGTs_recipient_genes.ffn'   % (pwd_MetaCHIP_op_folder, output_prefix)
        pwd_recipient_gene_seq_faa =    '%s/%s_x_detected_HGTs_recipient_genes.faa'   % (pwd_MetaCHIP_op_folder, output_prefix)
        pwd_donor_gene_seq_ffn =        '%s/%s_x_detected_HGTs_donor_genes.ffn'       % (pwd_MetaCHIP_op_folder, output_prefix)
        pwd_donor_gene_seq_faa =        '%s/%s_x_detected_HGTs_donor_genes.faa'       % (pwd_MetaCHIP_op_folder, output_prefix)

        pwd_detected_HGT_txt_handle = open(pwd_detected_HGT_txt, 'w')
        pwd_detected_HGT_txt_handle.write('Gene_1\tGene_2\tIdentity\tend_match\tfull_length_match\tdirection\n')
        recipient_gene_list = set()
        donor_gene_list = set()
        flanking_plot_file_list = set()
        for each_HGT in open(pwd_detected_HGT_PG_txt):

            if not each_HGT.startswith('Gene_1'):

                each_HGT_split = each_HGT.strip().split('\t')
                gene_1 = each_HGT_split[0]
                gene_2 = each_HGT_split[1]
                gene_1_genome = '_'.join(gene_1.split('_')[:-1])
                gene_2_genome = '_'.join(gene_2.split('_')[:-1])
                identity = float(each_HGT_split[4])
                end_match = each_HGT_split[5]
                full_length_match = each_HGT_split[6]
                direction = each_HGT_split[7]

                if direction != 'NA':
                    pwd_detected_HGT_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_1, gene_2, identity, end_match, full_length_match, direction))

                    recipient_genome = direction.split('-->')[1]
                    if gene_1_genome == recipient_genome:
                        recipient_gene_list.add(gene_1)
                        donor_gene_list.add(gene_2)
                    if gene_2_genome == recipient_genome:
                        recipient_gene_list.add(gene_2)
                        donor_gene_list.add(gene_1)

                    flanking_plot_file_list.add('%s___%s.SVG' % (gene_1, gene_2))

        pwd_detected_HGT_txt_handle.close()


        if len(recipient_gene_list) > 0:

            extract_donor_recipient_sequences(pwd_combined_ffn_file, recipient_gene_list, donor_gene_list,
                                              pwd_recipient_gene_seq_ffn, pwd_recipient_gene_seq_faa,
                                              pwd_donor_gene_seq_ffn, pwd_donor_gene_seq_faa)

            for each_flk_plot in flanking_plot_file_list:
                pwd_each_flk_plot = '%s/%s_x%s_Flanking_region_plots/1_Plots_normal/%s' % (pwd_MetaCHIP_op_folder, output_prefix, group_num, each_flk_plot)
                os.system('mv %s %s/%s_x%s_Flanking_region_plots/' % (pwd_each_flk_plot, pwd_MetaCHIP_op_folder, output_prefix, group_num))

            ###################################### Get_circlize_plot #######################################

            pwd_plot_circos =               '%s/%s_x%s_HGTs_among_provided_groups.pdf'                   % (pwd_MetaCHIP_op_folder, output_prefix, group_num)

            # get genome to group dict
            genome_to_group_dict = {}
            for genome in open(grouping_file):
                group_id2 = genome.strip().split(',')[0]
                genome_name = genome.strip().split(',')[1]
                genome_to_group_dict[genome_name] = group_id2

            Get_circlize_plot_customized_grouping(multi_level_detection, output_prefix, pwd_detected_HGT_txt, genome_to_group_dict, circos_HGT_R, pwd_plot_circos, pwd_MetaCHIP_op_folder)

            # remove tmp files
            os.remove(pwd_detected_HGT_PG_txt)
            os.system('rm -r %s/%s_x%s_Flanking_region_plots/1_Plots_normal'            % (pwd_MetaCHIP_op_folder, output_prefix, group_num))
            os.system('rm -r %s/%s_x%s_Flanking_region_plots/2_Plots_end_match'         % (pwd_MetaCHIP_op_folder, output_prefix, group_num))
            os.system('rm -r %s/%s_x%s_Flanking_region_plots/3_Plots_full_length_match' % (pwd_MetaCHIP_op_folder, output_prefix, group_num))
        else:
            report_and_log(('No HGT detected, program exited!'),pwd_log_file, keep_quiet)
            exit()

    else:
        detection_rank_list = args['r']

        # for single level detection
        if len(detection_rank_list) == 1:

            multi_level_detection = False
            pwd_MetaCHIP_op_folder_re = '%s/%s_HGT_ip%s_al%sbp_c%s_ei%s_f%skbp_%s*' % (MetaCHIP_wd, output_prefix, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp, detection_rank_list)

            MetaCHIP_op_folder = [os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)][0]
            group_num = int(MetaCHIP_op_folder.split('_')[-1][1:])

            pwd_MetaCHIP_op_folder =        '%s/%s'                                         % (MetaCHIP_wd, MetaCHIP_op_folder)
            pwd_detected_HGT_PG_txt =       '%s/%s_%s%s_HGTs_PG.txt'                        % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num)
            pwd_flanking_plot_folder =      '%s/%s_%s%s_Flanking_region_plots'              % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num)
            pwd_detected_HGT_txt =          '%s/%s_%s_detected_HGTs.txt'                    % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list)
            pwd_recipient_gene_seq_ffn =    '%s/%s_%s_detected_HGTs_recipient_genes.ffn'    % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list)
            pwd_recipient_gene_seq_faa =    '%s/%s_%s_detected_HGTs_recipient_genes.faa'    % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list)
            pwd_donor_gene_seq_ffn =        '%s/%s_%s_detected_HGTs_donor_genes.ffn'        % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list)
            pwd_donor_gene_seq_faa =        '%s/%s_%s_detected_HGTs_donor_genes.faa'        % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list)

            pwd_detected_HGT_txt_handle = open(pwd_detected_HGT_txt, 'w')
            pwd_detected_HGT_txt_handle.write('Gene_1\tGene_2\tIdentity\tend_match\tfull_length_match\tdirection\n')
            recipient_gene_list = set()
            donor_gene_list = set()
            flanking_plot_file_list = set()
            for each_HGT in open(pwd_detected_HGT_PG_txt):

                if not each_HGT.startswith('Gene_1'):

                    each_HGT_split = each_HGT.strip().split('\t')
                    gene_1 = each_HGT_split[0]
                    gene_2 = each_HGT_split[1]
                    gene_1_genome = '_'.join(gene_1.split('_')[:-1])
                    gene_2_genome = '_'.join(gene_2.split('_')[:-1])
                    identity = float(each_HGT_split[4])
                    end_match = each_HGT_split[5]
                    full_length_match = each_HGT_split[6]
                    direction = each_HGT_split[7]

                    if direction != 'NA':
                        pwd_detected_HGT_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_1, gene_2, identity, end_match, full_length_match, direction))

                        recipient_genome = direction.split('-->')[1]
                        if gene_1_genome == recipient_genome:
                            recipient_gene_list.add(gene_1)
                            donor_gene_list.add(gene_2)
                        if gene_2_genome == recipient_genome:
                            recipient_gene_list.add(gene_2)
                            donor_gene_list.add(gene_1)

                        flanking_plot_file_list.add('%s___%s.SVG' % (gene_1, gene_2))

            pwd_detected_HGT_txt_handle.close()

            extract_donor_recipient_sequences(pwd_combined_ffn_file, recipient_gene_list, donor_gene_list, pwd_recipient_gene_seq_ffn, pwd_recipient_gene_seq_faa, pwd_donor_gene_seq_ffn, pwd_donor_gene_seq_faa)

            for each_flk_plot in flanking_plot_file_list:
                pwd_each_flk_plot = '%s/%s_%s%s_Flanking_region_plots/1_Plots_normal/%s' % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num, each_flk_plot)
                os.system('mv %s %s/%s_%s%s_Flanking_region_plots/' % (pwd_each_flk_plot, pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num))

            ###################################### Get_circlize_plot #######################################

            grouping_file_re = '%s/%s_grouping_%s*.txt' % (MetaCHIP_wd, output_prefix, detection_rank_list)
            grouping_file = [os.path.basename(file_name) for file_name in glob.glob(grouping_file_re)][0]
            taxon_rank_num = grouping_file[len(output_prefix) + 1:].split('_')[0]
            pwd_grouping_file =         '%s/%s'                                     % (MetaCHIP_wd, grouping_file)
            pwd_plot_circos =           '%s/%s_%s_HGTs_among_provided_groups.pdf'   % (pwd_MetaCHIP_op_folder, output_prefix, taxon_rank_num)

            taxon_to_group_id_dict = {}
            for group in open(pwd_grouping_file):
                group_id = group.strip().split(',')[0]
                group_taxon = group.strip().split(',')[2]
                if group_id not in taxon_to_group_id_dict:
                    taxon_to_group_id_dict[group_id] = group_taxon

            # get genome to taxon dict
            genome_to_taxon_dict = {}
            for genome in open(pwd_grouping_file):
                group_id2 = genome.strip().split(',')[0]
                genome_name = genome.strip().split(',')[1]
                genome_to_taxon_dict[genome_name] = taxon_to_group_id_dict[group_id2]

            detected_HGT_num = -1
            for each_line in open(pwd_detected_HGT_txt):
                detected_HGT_num += 1

            if detected_HGT_num < 1:
                report_and_log(('No HGT detected, program exited!'), pwd_log_file, keep_quiet)
                exit()
            else:
                Get_circlize_plot(multi_level_detection, output_prefix, pwd_detected_HGT_txt, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, detection_rank_list, taxon_rank_num, pwd_MetaCHIP_op_folder)


            # remove tmp files
            os.remove(pwd_detected_HGT_PG_txt)
            os.system('rm -r %s/%s_%s%s_Flanking_region_plots/1_Plots_normal'               % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num))
            os.system('rm -r %s/%s_%s%s_Flanking_region_plots/2_Plots_end_match'            % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num))
            os.system('rm -r %s/%s_%s%s_Flanking_region_plots/3_Plots_full_length_match'    % (pwd_MetaCHIP_op_folder, output_prefix, detection_rank_list, group_num))


        # for multiple level detection
        if len(detection_rank_list) > 1:

            time_format = '[%Y-%m-%d %H:%M:%S]'
            print('%s Combine multiple level predictions' % (datetime.now().strftime(time_format)))

            multi_level_detection = True

            pwd_detected_HGT_txt_list = []
            pwd_flanking_plot_folder_list = []
            for detection_rank in detection_rank_list:
                if detection_rank not in ignored_rank_list:
                    pwd_MetaCHIP_op_folder_re = '%s/%s_HGT_ip%s_al%sbp_c%s_ei%s_f%skbp_%s*' % (MetaCHIP_wd, output_prefix, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp, detection_rank)
                    MetaCHIP_op_folder_list = [os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)]

                    if 'combined' not in MetaCHIP_op_folder_list[0]:
                        MetaCHIP_op_folder = MetaCHIP_op_folder_list[0]
                    else:
                        MetaCHIP_op_folder = MetaCHIP_op_folder_list[1]

                    group_num = int(MetaCHIP_op_folder.split('_')[-1][1:])
                    pwd_detected_HGT_txt        = '%s/%s/%s_%s%s_HGTs_PG.txt' % (MetaCHIP_wd, MetaCHIP_op_folder, output_prefix, detection_rank, group_num)
                    pwd_flanking_plot_folder    = '%s/%s/%s_%s%s_Flanking_region_plots' % (MetaCHIP_wd, MetaCHIP_op_folder, output_prefix, detection_rank, group_num)
                    pwd_detected_HGT_txt_list.append(pwd_detected_HGT_txt)
                    pwd_flanking_plot_folder_list.append(pwd_flanking_plot_folder)

            pwd_combined_prediction_folder = '%s/%s_combined_%s_HGTs_ip%s_al%sbp_c%s_ei%s_f%skbp' % (MetaCHIP_wd, output_prefix, detection_rank_list, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)

            genome_size_file =                          '%s/%s_all_genome_size.txt'                     % (MetaCHIP_wd, output_prefix)
            pwd_detected_HGT_txt_combined =             '%s/%s_%s_detected_HGTs.txt'                    % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
            pwd_recipient_gene_seq_ffn =                '%s/%s_%s_detected_HGTs_recipient_genes.ffn'    % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
            pwd_recipient_gene_seq_faa =                '%s/%s_%s_detected_HGTs_recipient_genes.faa'    % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
            pwd_donor_gene_seq_ffn =                    '%s/%s_%s_detected_HGTs_donor_genes.ffn'        % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
            pwd_donor_gene_seq_faa =                    '%s/%s_%s_detected_HGTs_donor_genes.faa'        % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
            pwd_flanking_plot_folder_combined_tmp =     '%s/%s_%s_Flanking_region_plots_tmp'            % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
            pwd_flanking_plot_folder_combined =         '%s/%s_%s_Flanking_region_plots'                % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)

            # combine prediction
            force_create_folder(pwd_combined_prediction_folder)
            combine_PG_output(pwd_detected_HGT_txt_list, output_prefix, detection_rank_list, pwd_detected_HGT_txt_combined)


            ############################################### extract sequences ##############################################

            # get recipient and donor gene list
            recipient_gene_list = set()
            recipient_genome_list = []
            donor_gene_list = set()
            plot_file_list = set()
            for each in open(pwd_detected_HGT_txt_combined):
                if not each.startswith('Gene_1'):

                    each_split = each.strip().split('\t')
                    gene_1 = each_split[0]
                    gene_2 = each_split[1]
                    gene_1_genome = '_'.join(gene_1.split('_')[:-1])
                    gene_2_genome = '_'.join(gene_2.split('_')[:-1])
                    direction = each_split[6]
                    plot_file = '%s___%s.SVG' % (gene_1, gene_2)
                    plot_file_list.add(plot_file)

                    recipient_genome = direction.split('-->')[1]
                    if '%)' in recipient_genome:
                        recipient_genome = recipient_genome.split('(')[0]
                    recipient_genome_list.append(recipient_genome)

                    if gene_1_genome == recipient_genome:
                        recipient_gene_list.add(gene_1)
                        donor_gene_list.add(gene_2)
                    if gene_2_genome == recipient_genome:
                        recipient_gene_list.add(gene_2)
                        donor_gene_list.add(gene_1)

            extract_donor_recipient_sequences(pwd_combined_ffn_file, recipient_gene_list, donor_gene_list, pwd_recipient_gene_seq_ffn, pwd_recipient_gene_seq_faa, pwd_donor_gene_seq_ffn, pwd_donor_gene_seq_faa)


            ############################################ combine flanking plots ############################################

            # create plot folders
            os.mkdir(pwd_flanking_plot_folder_combined_tmp)
            os.mkdir(pwd_flanking_plot_folder_combined)

            for flanking_plot_folder in pwd_flanking_plot_folder_list:
                flanking_plot_re = '%s/1_Plots_normal/*.SVG' % flanking_plot_folder
                flanking_plot_list = [os.path.basename(file_name) for file_name in glob.glob(flanking_plot_re)]

                for flanking_plot in flanking_plot_list:
                    pwd_flanking_plot = '%s/1_Plots_normal/%s' % (flanking_plot_folder, flanking_plot)
                    os.system('cp %s %s/' % (pwd_flanking_plot, pwd_flanking_plot_folder_combined_tmp))

                # os.system('cp %s/1_Plots_normal/* %s/' % (flanking_plot_folder, pwd_flanking_plot_folder_combined_tmp))

            for plot_file in plot_file_list:
                pwd_plot_file = '%s/%s' % (pwd_flanking_plot_folder_combined_tmp, plot_file)
                os.system('mv %s %s/' % (pwd_plot_file, pwd_flanking_plot_folder_combined))


            ###################################### Get_circlize_plot #######################################

            for detection_rank in detection_rank_list:
                if detection_rank not in ignored_rank_list:
                    grouping_file_re = '%s/%s_grouping_%s*.txt' % (MetaCHIP_wd, output_prefix, detection_rank)
                    grouping_file = [os.path.basename(file_name) for file_name in glob.glob(grouping_file_re)][0]
                    taxon_rank_num = grouping_file.split('.')[-2].split('_')[-1][1:]
                    pwd_grouping_file =         '%s/%s'                         % (MetaCHIP_wd, grouping_file)
                    pwd_plot_circos =           '%s/%s_HGTs_among_%s.pdf'                   % (pwd_combined_prediction_folder, output_prefix, rank_abbre_dict_plural[detection_rank])
                    pwd_df_circos =             '%s/%s_HGTs_among_%s.tab'                   % (pwd_combined_prediction_folder, output_prefix, rank_abbre_dict_plural[detection_rank])

                    # get genome to taxon dict
                    genome_to_taxon_dict = {}
                    for genome in open(pwd_grouping_file):
                        genome_name = genome.strip().split(',')[1]
                        genome_taxon = genome.strip().split(',')[2]
                        genome_to_taxon_dict[genome_name] = genome_taxon

                    detected_HGT_num = -1
                    for each_line in open(pwd_detected_HGT_txt_combined):
                        detected_HGT_num += 1

                    if detected_HGT_num < 1:
                        report_and_log(('No HGT detected, program exited!'), pwd_log_file, keep_quiet)
                        exit()
                    else:
                        Get_circlize_plot(multi_level_detection, output_prefix, pwd_detected_HGT_txt_combined, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, detection_rank, taxon_rank_num, pwd_combined_prediction_folder)

            ###################################### remove tmp files #######################################

            # remove tmp files
            os.system('rm -r %s' % pwd_flanking_plot_folder_combined_tmp)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-o',             required=False, default=None,                 help='output folder (default: current working directory)')
    parser.add_argument('-p',             required=True,                                help='output prefix')
    parser.add_argument('-r',             required=False, default=None,                 help='grouping rank')
    parser.add_argument('-g',             required=False, default=None,                 help='grouping file')
    parser.add_argument('-cov',           required=False, type=int,     default=75,     help='coverage cutoff, default: 75')
    parser.add_argument('-al',            required=False, type=int,     default=200,    help='alignment length cutoff, default: 200')
    parser.add_argument('-flk',           required=False, type=int,     default=10,     help='the length of flanking sequences to plot (Kbp), default: 10')
    parser.add_argument('-ip',            required=False, type=int,     default=90,     help='identity percentile cutoff, default: 90')
    parser.add_argument('-ei',            required=False, type=float,   default=80,     help='end match identity cutoff, default: 80')
    parser.add_argument('-t',             required=False, type=int,     default=1,      help='number of threads, default: 1')
    parser.add_argument('-NoEbCheck',     required=False, action="store_true",          help='disable end break and contig match check for fast processing, not recommend for metagenome-assembled genomes (MAGs)')
    parser.add_argument('-force',         required=False, action="store_true",          help='overwrite previous results')
    parser.add_argument('-quiet',         required=False, action="store_true",          help='Do not report progress')
    parser.add_argument('-tmp',           required=False, action="store_true",          help='keep temporary files')

    args = vars(parser.parse_args())

    BP(args, config_dict)

