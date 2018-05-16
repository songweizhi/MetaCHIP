#!/usr/bin/env python
import os
import shutil
import argparse
import itertools
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
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from string import ascii_uppercase
from datetime import datetime


# requirement for contig id
# bin name splitter for direction information
# check the m0 group for the big differentces
# if run_blast = 1, there is no need to provide the blast results
# endbreak checked twice, remove one


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
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Filtering blast matches with the following criteria: Query name != Subject name, Alignment length >= %sbp and coverage >= %s%s' % (alignment_length_cutoff, coverage_cutoff, '%'))
    for match in matches:
        #print('\r%s x 1000 blast matches detected in total, filtering the %s x 1000th' % (float("{0:.2f}".format(total_match_number/1000)), float("{0:.2f}".format(n/1000))))
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
            if (query_bin_name != subject_bin_name) and (align_len >= int(alignment_length_cutoff)) and (coverage_q >= int(coverage_cutoff)) and (coverage_s >= int(coverage_cutoff)):
                out_temp.write(match)
                # remove the same match but with swapped query-subject position
                if (query_name_subject_name not in counted_match) and (subject_name_query_name not in counted_match):
                    all_identities.append(identity)
                    counted_match.append(query_name_subject_name)
        n += 1

    out_temp.close()
    return all_identities


def do():
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

    # check length
    if len(current_group_pair_identities) >= minimum_plot_number:
        if current_group_pair_name == current_group_pair_name_swapped:
            plot_identity_list(current_group_pair_identities, 'None', current_group_pair_name, pwd_iden_distrib_plot_folder)
        else:
            plot_identity_list(current_group_pair_identities, current_group_pair_identity_cut_off, current_group_pair_name, pwd_iden_distrib_plot_folder)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " Plotting identity distribution (%dth): %s" % (ploted_group, current_group_pair_name,))
    else:
        unploted_groups.write('%s\t%s\n' % (current_group_pair_name, len(current_group_pair_identities)))
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " Plotting identity distribution (%dth): %s, blast hits < %d, skipped" % (ploted_group, current_group_pair_name, minimum_plot_number))


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


def plot_identity_list(identity_list, identity_cut_off, title, output_foler) :
    identity_list = sorted(identity_list)
    # get statistics
    match_number = len(identity_list)
    average_iden = float(np.average(identity_list))
    average_iden = float("{0:.2f}".format(average_iden))
    max_match = float(np.max(identity_list))
    min_match = float(np.min(identity_list))
    # get hist plot
    num_bins = 50
    plt.hist(identity_list,
             num_bins,
             alpha = 0.1,
             normed = 1,
             facecolor = 'blue')  # normed = 1 normalized to 1, that is probablity
    plt.title('Group: %s' % title)
    plt.xlabel('Identity')
    plt.ylabel('Probability')
    plt.subplots_adjust(left = 0.15)
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
                xy = (identity_cut_off, 0),
                xytext = (identity_cut_off, density(identity_cut_off)),
                arrowprops = dict(width = 0.5,
                                  headwidth = 0.5,
                                  facecolor = 'red',
                                  edgecolor = 'red',
                                  shrink = 0.02))
    else:
        pass
    # Get plot
    plt.savefig('%s/%s.png' % (output_foler, title), dpi = 300)
    plt.close()


def get_candidates(targets_group_file, gene_with_g_file_name, gene_only_name_file_name, group_pair_iden_cutoff_dict):
    targets_group = open(targets_group_file)
    output_1 = open(gene_with_g_file_name, 'w')
    output_2 = open(gene_only_name_file_name, 'w')

    for group in targets_group:
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


def set_contig_track_features(gene_contig, name_group_dict, candidate_list, feature_set):
    # add features to feature set
    for feature in gene_contig.features:
        if feature.type == "CDS":
            # define label color
            if feature.qualifiers['locus_tag'][0] in candidate_list:
                label_color = colors.blue
            else:
                label_color = colors.black

            # change gene name
            bin_name_gbk_split = feature.qualifiers['locus_tag'][0].split('_')
            bin_name_gbk = '_'.join(bin_name_gbk_split[:-1])
            feature.qualifiers['locus_tag'][0] = name_group_dict[bin_name_gbk] + '_' + bin_name_gbk_split[-1]

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
                                    arrowshaft_height=0.7,
                                    arrowhead_length=0.4,
                                    label_color=label_color,
                                    label_size=10,
                                    label_angle=label_angle,
                                    label_position="middle")


def get_flanking_region(input_gbk_file, HGT_candidate, flanking_length):

    wd, gbk_file = os.path.split(input_gbk_file)
    new_gbk_file = '%s/%s_%sbp_temp.gbk' % (wd, HGT_candidate, flanking_length)
    new_gbk_final_file = '%s/%s_%sbp.gbk' % (wd, HGT_candidate, flanking_length)
    new_fasta_final_file = '%s/%s_%sbp.fasta' % (wd, HGT_candidate, flanking_length)
    output_plot = '%s/%s_%sbp.eps' % (wd, HGT_candidate, flanking_length)

    # get flanking range of candidate
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_start = 0
    new_end = 0
    contig_length = 0
    for record in input_gbk:
        for gene in record.features:
            # get contig length
            if gene.type == 'source':
                contig_length = int(gene.location.end)
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
                    #print(record.seq[gene.location.start:gene.location.end])
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
                               id = record.id,
                               name = record.name,
                               description = record.description,
                               annotations = record.annotations)

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
                #print(record.seq[gene.location.start:gene.location.end])
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


def get_end_break_value(dict_value_list_in, ending_match_length):
    gene_1_start = dict_value_list_in[0][1]
    gene_1_end = dict_value_list_in[0][2]
    gene_1_ctg_length = dict_value_list_in[0][4]

    gene_2_start = dict_value_list_in[1][1]
    gene_2_end = dict_value_list_in[1][2]
    gene_2_ctg_length = dict_value_list_in[1][4]

    # for gene1
    gene1_end_location = 0
    if (gene_1_start <= ending_match_length) or ((gene_1_ctg_length - gene_1_end) <= ending_match_length):
        gene1_end_location = 1

    # for gene2
    gene2_end_location = 0
    if (gene_2_start <= ending_match_length) or ((gene_2_ctg_length - gene_2_end) <= ending_match_length):
        gene2_end_location = 1

    if (gene1_end_location == 1) and (gene2_end_location == 1):
        end_break = 1
    else:
        end_break = 0

    return end_break


def check_end_break(folder_name, flanking_length, end_seq_length, pwd_blastn_exe):

    # define file name
    recipient_gene = folder_name.split('___')[0]
    donor_gene = folder_name.split('___')[1]
    file_recipient_gene_3000_gbk = '%s_%sbp.gbk' % (recipient_gene, flanking_length)
    file_donor_gene_3000_gbk = '%s_%sbp.gbk' % (donor_gene, flanking_length)

    # read in recipient/donor contig
    recipient_contig_record = SeqIO.read(file_recipient_gene_3000_gbk, 'genbank')
    recipient_contig_seq = recipient_contig_record.seq
    donor_contig_record = SeqIO.read(file_donor_gene_3000_gbk, 'genbank')
    donor_contig_seq = donor_contig_record.seq

    # get ending sequence of the recipient and donor genes
    ending_seq_description = ''

    # export recipient_left_end_seq
    recipient_left_end_seq = recipient_contig_seq[0:end_seq_length]
    recipient_left_end_id = '%s_le%s' % (recipient_gene, end_seq_length)
    recipient_left_end_handle = '%s/%s.fasta' % (os.getcwd(), recipient_left_end_id)
    export_dna_record(recipient_left_end_seq, recipient_left_end_id, ending_seq_description, recipient_left_end_handle)

    # export recipient_right_end_seq
    recipient_right_end_seq = recipient_contig_seq[len(recipient_contig_seq) - end_seq_length:]
    recipient_right_end_id = '%s_re%s' % (recipient_gene, end_seq_length)
    recipient_right_end_handle = '%s/%s.fasta' % (os.getcwd(), recipient_right_end_id)
    export_dna_record(recipient_right_end_seq, recipient_right_end_id, ending_seq_description, recipient_right_end_handle)

    # export donor_left_end_seq
    donor_left_end_seq = donor_contig_seq[0:end_seq_length]
    donor_left_end_id = '%s_le%s' % (donor_gene, end_seq_length)
    donor_left_end_handle = '%s/%s.fasta' % (os.getcwd(), donor_left_end_id)
    export_dna_record(donor_left_end_seq, donor_left_end_id, ending_seq_description, donor_left_end_handle)

    # export donor_right_end_seq
    donor_right_end_seq = donor_contig_seq[len(donor_contig_seq) - end_seq_length:]
    donor_right_end_id = '%s_re%s' % (donor_gene, end_seq_length)
    donor_right_end_handle = '%s/%s.fasta' % (os.getcwd(), donor_right_end_id)
    export_dna_record(donor_right_end_seq, donor_right_end_id, ending_seq_description, donor_right_end_handle)

    # run blastn between ending sequences:
    blast_parameters = '-evalue 1e-5 -outfmt 6 -task blastn'
    output_rle_dle = '%s/%s_rle___%s_dle.tab' % (os.getcwd(), recipient_gene, donor_gene)
    output_rle_dre = '%s/%s_rle___%s_dre.tab' % (os.getcwd(), recipient_gene, donor_gene)
    output_rre_dle = '%s/%s_rre___%s_dle.tab' % (os.getcwd(), recipient_gene, donor_gene)
    output_rre_dre = '%s/%s_rre___%s_dre.tab' % (os.getcwd(), recipient_gene, donor_gene)
    compare_end_blast_rle_dle = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_left_end_handle, donor_left_end_handle, output_rle_dle, blast_parameters)
    compare_end_blast_rle_dre = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_left_end_handle, donor_right_end_handle, output_rle_dre, blast_parameters)
    compare_end_blast_rre_dle = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_right_end_handle, donor_left_end_handle, output_rre_dle, blast_parameters)
    compare_end_blast_rre_dre = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_right_end_handle, donor_right_end_handle, output_rre_dre, blast_parameters)
    os.system(compare_end_blast_rle_dle)
    os.system(compare_end_blast_rle_dre)
    os.system(compare_end_blast_rre_dle)
    os.system(compare_end_blast_rre_dre)

    match_files = [output_rle_dle, output_rle_dre, output_rre_dle, output_rre_dre]
    match_profile = []
    for each_match in match_files:
        match_dir_list = []
        for each_hit in open(each_match):
            each_hit_split = each_hit.strip().split('\t')
            qstart = int(each_hit_split[6])
            qend = int(each_hit_split[7])
            sstart = int(each_hit_split[8])
            send = int(each_hit_split[9])
            q_direction = qend - qstart
            s_direction = send - sstart
            q_dir = ''
            if q_direction > 0:
                q_dir = 'forward'
            if q_direction < 0:
                q_dir = 'backward'
            s_dir = ''
            if s_direction > 0:
                s_dir = 'forward'
            if s_direction < 0:
                s_dir = 'backward'
            match_dir = ''
            if q_dir == s_dir:
                match_dir = 'same direction'
            if q_dir != s_dir:
                match_dir = 'opposite direction'
            match_dir_list.append(match_dir)
        match_dir_list_uniq = []
        for each in match_dir_list:
            if each not in match_dir_list_uniq:
                match_dir_list_uniq.append(each)
        match_profile.append(match_dir_list_uniq)

    break_end = ''
    if (match_profile[0] == ['opposite direction']) or (match_profile[1] == ['same direction']) or (match_profile[2] == ['same direction']) or (match_profile[3] == ['opposite direction']):
        break_end = True
    else:
        break_end = False

    return break_end


def get_gbk_blast_act(candidates_file, gbk_file, flanking_length, end_seq_length, name_to_group_number_dict, path_to_output_act_folder, pwd_blastn_exe, keep_temp):

    matches = open(candidates_file)

    total = 0
    for match in matches:
        total += 1

    n = 1
    candidates_2_endbreak_dict = {}
    candidates_2_endlocation_dict = {}

    for match in open(candidates_file):
        genes = match.strip().split('\t')[:-1]
        folder_name = '___'.join(genes)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " Plotting flanking region (%d/%d): %s" % (n, total, folder_name))
        os.mkdir('%s/%s' % (path_to_output_act_folder, folder_name))

        dict_value_list = []
        # Extract gbk and faa files
        records = SeqIO.parse(gbk_file, 'genbank')
        for record in records:
            for gene_f in record.features:
                if 'locus_tag' in gene_f.qualifiers:

                    for gene_1 in genes:
                        if gene_1 in gene_f.qualifiers["locus_tag"]:
                            dict_value_list.append([gene_1, int(gene_f.location.start), int(gene_f.location.end), gene_f.location.strand, len(record.seq)])
                            pwd_gbk_file = '%s/%s/%s.gbk' % (path_to_output_act_folder, folder_name, gene_1)
                            pwd_fasta_file = '%s/%s/%s.fasta' % (path_to_output_act_folder, folder_name, gene_1)
                            SeqIO.write(record, pwd_gbk_file, 'genbank')
                            SeqIO.write(record, pwd_fasta_file, 'fasta')
                            # get flanking regions
                            get_flanking_region(pwd_gbk_file, gene_1, flanking_length)

        candidates_2_endlocation_dict[folder_name] = dict_value_list

        # Run Blast
        prefix_c = '%s/%s' % (path_to_output_act_folder, folder_name)
        output_c = '%s/%s.txt' % (prefix_c, folder_name)
        query_c = '%s/%s_%sbp.fasta' % (prefix_c, genes[0], flanking_length)
        subject_c = '%s/%s_%sbp.fasta' % (prefix_c, genes[1], flanking_length)

        parameters_c_n = '-evalue 1e-5 -outfmt 6 -task blastn'
        command_blast = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_c, subject_c, output_c, parameters_c_n)
        os.system(command_blast)

        # read in gbk files
        matche_pair_list = []
        for each_gene in genes:
            path_to_gbk_file = '%s/%s/%s_%sbp.gbk' % (path_to_output_act_folder, folder_name, each_gene, flanking_length)
            gene_contig = SeqIO.read(path_to_gbk_file, "genbank")
            matche_pair_list.append(gene_contig)
        bin_record_list = []
        bin_record_list.append(matche_pair_list)

        # create an empty diagram
        diagram = GenomeDiagram.Diagram()

        # add tracks to diagram
        max_len = 0
        for gene1_contig, gene2_contig in bin_record_list:
            # set diagram track length
            max_len = max(max_len, len(gene1_contig), len(gene2_contig))

            # add gene content track for gene1_contig
            contig_1_gene_content_track = diagram.new_track(1,
                                                            name=gene1_contig.name,
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
                                                            name=gene2_contig.name,
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
            set_contig_track_features(gene1_contig, name_to_group_number_dict, genes, feature_sets_1)
            set_contig_track_features(gene2_contig, name_to_group_number_dict, genes, feature_sets_2)

            ####################################### add crosslink from blast results #######################################

            path_to_blast_result = '%s/%s/%s.txt' % (path_to_output_act_folder, folder_name, folder_name)
            blast_results = open(path_to_blast_result)

            # parse blast results
            for each_line in blast_results:
                each_line_split = each_line.split('\t')
                query = each_line_split[0]
                identity = float(each_line_split[2])
                query_start = int(each_line_split[6])
                query_end = int(each_line_split[7])
                target_start = int(each_line_split[8])
                target_end = int(each_line_split[9])

                # use color to reflect identity
                color = colors.linearlyInterpolatedColor(colors.white, colors.red, 50, 100, identity)

                # determine which is which (query/target to contig_1/contig_2)
                # if query is contig_1
                if query == gene1_contig.name :
                    link = CrossLink((contig_1_gene_content_track, query_start, query_end),
                                     (contig_2_gene_content_track, target_start, target_end),
                                     color = color,
                                     border = color,
                                     flip = False)
                    diagram.cross_track_links.append(link)

                # if query is contig_2
                elif query == gene2_contig.name:
                    link = CrossLink((contig_2_gene_content_track, query_start, query_end),
                                     (contig_1_gene_content_track, target_start, target_end),
                                     color = color,
                                     border = color,
                                     flip = False)
                    diagram.cross_track_links.append(link)

            ############################################### Draw and Export ################################################

            diagram.draw(format="linear",
                         orientation="landscape",
                         pagesize=(75 * cm, 25 * cm),
                         fragments=1,
                         start=0,
                         end=max_len)

            diagram.write('%s/%s.eps' % (path_to_output_act_folder, folder_name), "EPS")

        # get match category
        current_wd = os.getcwd()
        os.chdir('%s/%s' % (path_to_output_act_folder, folder_name))

        end_break = get_end_break_value(dict_value_list, end_seq_length)
        candidates_2_endbreak_dict[folder_name] = end_break
        os.chdir(path_to_output_act_folder)

        if end_break == True:
            os.system('mv %s.eps 0_end_break/' % folder_name)
            os.chdir(current_wd)

        # remove temporary folder
        os.chdir(current_wd)
        if keep_temp == 0:
            shutil.rmtree('%s/%s' % (path_to_output_act_folder, folder_name), ignore_errors=True)
            if os.path.isdir('%s/%s' % (path_to_output_act_folder, folder_name)):
                shutil.rmtree('%s/%s' % (path_to_output_act_folder, folder_name), ignore_errors=True)
        n += 1

    return candidates_2_endbreak_dict


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


parser = argparse.ArgumentParser()

parser.add_argument('-g',
                    required=True,
                    help='grouping file')

parser.add_argument('-ffn',
                    required=True,
                    help='combined ffn file')

parser.add_argument('-gbk',
                    required=True,
                    help='combined gbk file')

parser.add_argument('-blastall',
                    required=True,
                    help='all vs all blast results')

parser.add_argument('-cov',
                    required=False,
                    type=int,
                    default=70,
                    help='coverage cutoff')

parser.add_argument('-al',
                    required=False,
                    type=int,
                    default=200,
                    help='alignment length cutoff')

parser.add_argument('-flk',
                    required=False,
                    type=int,
                    default=3000,
                    help='the length of flanking sequences to plot')

parser.add_argument('-ip',
                    required=False,
                    type=int,
                    default=90,
                    help='identity percentile')

parser.add_argument('-eb',
                    required=False,
                    type=int,
                    default=1000,
                    help='the minimal length to be considered as end break')

parser.add_argument('-tmp',
                    action="store_true",
                    required=False,
                    help='keep temporary files')

parser.add_argument('-blast',
                    action="store_true",
                    required=False,
                    help='run blast step or not')

parser.add_argument('-blastn',
                    required=False,
                    default='blastn',
                    help='path to blastn executable')

parser.add_argument('-makeblastdb',
                    required=False,
                    default='makeblastdb',
                    help='path to makeblastdb executable')

args = vars(parser.parse_args())

grouping_file = args['g']
ffn_file = args['ffn']
gbk_file = args['gbk']
blast_results = args['blastall']
cover_cutoff = args['cov']
flanking_length = args['flk']
identity_percentile = args['ip']
align_len_cutoff = args['al']
ending_match_length = args['eb']
run_blastn = args['blast']
keep_temp = args['tmp']
pwd_blastn_exe = args['blastn']
pwd_makeblastdb_exe = args['makeblastdb']


############################################### Define folder/file name ################################################

print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Define folder/file names and create output folder')

wd = os.getcwd()
pwd_grouping_file = '%s/%s' % (wd, grouping_file)
op_folder = 'output_ip%s_al%sbp_c%s_e%sbp_g%s' % (str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(ending_match_length), get_number_of_group(pwd_grouping_file))

iden_distrib_plot_folder =                          'identity_distribution'
qual_idens_file =                                   'qualified_identities.txt'
qual_idens_file_gg =                                'qualified_identities_gg.txt'
qual_idens_file_gg_sorted =                         'qualified_identities_gg_sorted.txt'
unploted_groups_file =                              'unploted_groups.txt'
qual_idens_with_group_filename =                    'qualified_identities_with_group.txt'
qual_idens_with_group_sorted_filename =             'qualified_identities_with_group_sorted.txt'
subjects_in_one_line_filename =                     'subjects_in_one_line.txt'
group_pair_iden_cutoff_file_name =                  'identity_cutoff.txt'
op_candidates_with_group_file_name =                'HGT_candidates_with_group.txt'
op_candidates_only_gene_file_name =                 'HGT_candidates_only_id.txt'
op_candidates_only_gene_file_name_uniq =            'HGT_candidates_uniq.txt'
op_candidates_only_gene_uniq_end_break =            'HGT_candidates.txt'
op_candidates_seq_nc =                              'HGT_candidates_nc.fasta'
op_candidates_seq_aa =                              'HGT_candidates_aa.fasta'
op_act_folder_name =                                'Flanking_regions'
gbk_subset_file =                                   'combined_subset.gbk'

pwd_iden_distrib_plot_folder =                      '%s/%s/%s'    % (wd, op_folder, iden_distrib_plot_folder)
pwd_qual_iden_file =                                '%s/%s/%s'    % (wd, op_folder, qual_idens_file)
pwd_qual_iden_file_gg =                             '%s/%s/%s'    % (wd, op_folder, qual_idens_file_gg)
pwd_qual_iden_file_gg_sorted =                      '%s/%s/%s'    % (wd, op_folder, qual_idens_file_gg_sorted)
pwd_unploted_groups_file =                          '%s/%s/%s/%s' % (wd, op_folder, iden_distrib_plot_folder, unploted_groups_file)
pwd_qual_idens_with_group =                         '%s/%s/%s'    % (wd, op_folder, qual_idens_with_group_filename)
pwd_qual_idens_with_group_sorted =                  '%s/%s/%s'    % (wd, op_folder, qual_idens_with_group_sorted_filename)
pwd_subjects_in_one_line =                          '%s/%s/%s'    % (wd, op_folder, subjects_in_one_line_filename)
pwd_group_pair_iden_cutoff_file =                   '%s/%s/%s'    % (wd, op_folder, group_pair_iden_cutoff_file_name)
pwd_op_candidates_with_group_file =                 '%s/%s/%s'    % (wd, op_folder, op_candidates_with_group_file_name)
pwd_op_candidates_only_gene_file =                  '%s/%s/%s'    % (wd, op_folder, op_candidates_only_gene_file_name)
pwd_op_candidates_only_gene_file_uniq =             '%s/%s/%s'    % (wd, op_folder, op_candidates_only_gene_file_name_uniq)
pwd_op_cans_only_gene_uniq_end_break =              '%s/%s/%s'    % (wd, op_folder, op_candidates_only_gene_uniq_end_break)
pwd_op_candidates_seq_nc =                          '%s/%s/%s'    % (wd, op_folder, op_candidates_seq_nc)
pwd_op_candidates_seq_aa =                          '%s/%s/%s'    % (wd, op_folder, op_candidates_seq_aa)
pwd_gbk_subset_file =                               '%s/%s/%s'    % (wd, op_folder, gbk_subset_file)
pwd_op_act_folder =                                 '%s/%s/%s'    % (wd, op_folder, op_act_folder_name)
pwd_grouping_file_with_id =                         '%s/%s/%s'    % (wd, op_folder, 'grouping_with_id.txt')
pwd_ffn_file =                                      '%s/%s'       % (wd, ffn_file)
pwd_gbk_file =                                      '%s/%s'       % (wd, gbk_file)
pwd_blast_results = ''
if run_blastn == 0:
    pwd_blast_results =                             '%s/%s'       % (wd, blast_results)
if run_blastn == 1:
    pwd_blast_results =                             '%s/%s'       % (wd, 'all_vs_all_ffn.tab')

# check whether file exist
# unfound_inputs = []
# for each_input in [pwd_grouping_file, pwd_prokka_output]:
#     if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
#         unfound_inputs.append(each_input)
#
# if run_blastn == 0:
#     if not os.path.isfile(pwd_blast_results):
#         unfound_inputs.append(pwd_blast_results)
#
# if len(unfound_inputs) > 0:
#     for each_unfound in unfound_inputs:
#         print('%s not found' % each_unfound)
#     exit()

########################################################################################################################

# # Prepare input files
# os.system('cat %s/*/*.gbk > %s' % (pwd_prokka_output, pwd_gbk_file)) # get combined gbk file
# os.system('cat %s/%s/*/*.ffn > combined.ffn' % (wd, prokka_output)) # get combined ffn file

# run blastn if specified
if run_blastn == 1:
    if os.path.isdir('blastdb'):
        shutil.rmtree('blastdb')
    os.mkdir('blastdb')

    os.system('cp %s ./blastdb/' % pwd_ffn_file)
    makeblastdb_cmd = '%s -in blastdb/%s -dbtype nucl -parse_seqids' % (ffn_file, pwd_makeblastdb_exe)
    os.system(makeblastdb_cmd)
    blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running blastn, be patient...')
    os.system('%s -query %s -db blastdb/%s -out %s %s' % (pwd_blastn_exe, pwd_ffn_file, ffn_file, pwd_blast_results, blast_parameters))

# create outputs folder
if os.path.isdir(op_folder):
    shutil.rmtree(op_folder, ignore_errors=True)
    if os.path.isdir(op_folder):
        shutil.rmtree(op_folder, ignore_errors=True)
    os.makedirs(op_folder)
else:
    os.makedirs(op_folder)


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

####################################################### Main code ######################################################

sleep(1)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Plotting identity distribution between each pair of groups')

# get qualified identities after alignment length and coverage filter (self-match excluded)
all_identities = get_all_identity_list(pwd_blast_results, genome_name_list, align_len_cutoff, cover_cutoff, pwd_qual_iden_file)

# create folder to hold group-group identity distribution plot
os.makedirs(pwd_iden_distrib_plot_folder)

# plot overall identity distribution
all_identities_plot_title = 'All_vs_All'
plot_identity_list(all_identities, 'None', all_identities_plot_title, pwd_iden_distrib_plot_folder)


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
            do()
            ploted_group += 1
            # restore current_group_pair_name and current_group_pair_identities for next group pair
            current_group_pair_name = group_pair
            current_group_pair_identities = []
            current_group_pair_identities.append(identity)
# for the last group
do()
group_pair_iden_cutoff_file.close()

sleep(1)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Analyzing Blast matches to get HGT candidates')

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

#################################################### Get ACT images ####################################################

sleep(1)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Preparing subset of combined.gbk file for flanking region plotting')

# get gene list of all candidates
all_candidates_genes = []
for match in open(pwd_op_candidates_only_gene_file):
    match_split = match.strip().split('\t')
    for gene in match_split:
        if gene not in all_candidates_genes:
            all_candidates_genes.append(gene)

# get subset of combined gbk file
gbk_subset = open(pwd_gbk_subset_file, 'w')
records_recorded = []
for record in SeqIO.parse(pwd_gbk_file, 'genbank'):
    record_id = record.id
    for gene_f in record.features:
        if 'locus_tag' in gene_f.qualifiers:
            for gene_r in all_candidates_genes:
                if gene_r in gene_f.qualifiers["locus_tag"]:
                    if record_id not in records_recorded:
                        SeqIO.write(record, gbk_subset, 'genbank')
                        records_recorded.append(record_id)
gbk_subset.close()


# create folder to hold ACT output
os.makedirs(pwd_op_act_folder)
os.makedirs('%s/0_end_break' % pwd_op_act_folder)

# plot flanking regions
candidates_2_endbreak_dict = get_gbk_blast_act(pwd_op_candidates_only_gene_file_uniq, pwd_gbk_subset_file, flanking_length, ending_match_length, name_to_group_number_dict, pwd_op_act_folder, pwd_blastn_exe, keep_temp)

# add end break information to output file
output_file = open(pwd_op_cans_only_gene_uniq_end_break, 'w')
output_file.write('Gene_1\tGene_2\tGene_1_group\tGene_2_group\tIdentity(nc)\tEnd_break\n')
all_candidate_genes = []
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
    end_break = candidates_2_endbreak_dict[concatenated]
    # write to output files
    if end_break == 1:
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_group, donor_genome_group, identity, 'yes' ))
    else:
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_group, donor_genome_group, identity, 'no'))
output_file.close()

# get qualified HGT candidates
qualified_HGT_candidates = []
for each_candidate_2 in open(pwd_op_cans_only_gene_uniq_end_break):
    each_candidate_2_split = each_candidate_2.strip().split('\t')
    end_break_2 = each_candidate_2_split[5]
    if end_break_2 == 'no':
        if each_candidate_2_split[0] not in qualified_HGT_candidates:
            qualified_HGT_candidates.append(each_candidate_2_split[0])
        if each_candidate_2_split[1] not in qualified_HGT_candidates:
            qualified_HGT_candidates.append(each_candidate_2_split[1])

# export nc and aa sequence of candidate HGTs
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Extracting nc and aa sequences for predicted HGTs')
candidates_seq_nc_handle = open(pwd_op_candidates_seq_nc, 'w')
candidates_seq_aa_handle = open(pwd_op_candidates_seq_aa, 'w')
for each_seq in SeqIO.parse(pwd_ffn_file, 'fasta'):
    if each_seq.id in qualified_HGT_candidates:
        SeqIO.write(each_seq, candidates_seq_nc_handle, 'fasta')
        each_seq_aa = each_seq
        each_seq_aa.seq = each_seq.seq.translate()
        SeqIO.write(each_seq_aa, candidates_seq_aa_handle, 'fasta')
candidates_seq_nc_handle.close()
candidates_seq_aa_handle.close()

# remove temporary files
if keep_temp == 0:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Deleting temporary files')
    os.remove(pwd_qual_iden_file)
    os.remove(pwd_qual_iden_file_gg)
    os.remove(pwd_qual_iden_file_gg_sorted)
    os.remove(pwd_qual_idens_with_group)
    os.remove(pwd_qual_idens_with_group_sorted)
    os.remove(pwd_subjects_in_one_line)
    os.remove(pwd_gbk_subset_file)
    os.remove(pwd_op_candidates_only_gene_file_uniq)
    os.remove(pwd_op_candidates_with_group_file)
    os.remove(pwd_op_candidates_only_gene_file)

sleep(1)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Done for Best-match approach prediction!')
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' You may want to run the Explicit tree approach for further validation.')
