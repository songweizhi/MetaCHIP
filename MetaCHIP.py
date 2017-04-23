#!/usr/bin/env python
import os
import shutil
import argparse
import configparser
from sys import stdout
from time import sleep
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm


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
    print('Filtering blast matches with the following criteria: Query name != Subject name, Alignment length >= %sbp and coverage >= %s%s' % (alignment_length_cutoff, coverage_cutoff, '%'))
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
            if (query_bin_name != subject_bin_name) and (align_len >= int(alignment_length_cutoff)) and (coverage_q >= int(coverage_cutoff)) and (coverage_s >= int(coverage_cutoff)):
                out_temp.write(match)
                # remove the same match but with swapped query-subject position
                if (query_name_subject_name in counted_match) or (subject_name_query_name in counted_match):
                    pass
                else:
                    all_identities.append(identity)
                    counted_match.append(query_name_subject_name)
            else:
                pass
        else:
            pass
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
            for each_subject in subjects_list:
                each_subject_g = each_subject[0]
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
                    if each_g[0] not in non_self_group_subject_list_uniq:
                        non_self_group_subject_list_uniq.append(each_g[0])
                    else:
                        pass

                # if all non-self-group subjects come from the same group
                if len(non_self_group_subject_list_uniq) == 1:

                    # get the maximum and average identity from self-group
                    sg_maximum = 0
                    sg_sum = 0
                    sg_subject_number = 0
                    for each_sg_subject in self_group_subject_list:
                        each_sg_subject_iden = float(each_sg_subject.split('|')[2])
                        if each_sg_subject_iden > sg_maximum :
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
                    else:
                        pass

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
                            each_nsg_subject_iden = float(each_nsg_subject.split('|')[2])
                            if each_nsg_subject[0] == each_nsg:
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
                        else:
                            pass

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
                        else:
                            pass
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
            if feature.location.strand == 1:
                label_angle = 45
                color = colors.lightblue
            elif feature.location.strand == -1:
                label_angle = -225
                color = colors.lightgreen
            # add feature
            feature_set.add_feature(feature,
                                    color = color,
                                    label = True,
                                    sigil = 'ARROW',
                                    arrowshaft_height = 0.7,
                                    arrowhead_length = 0.4,
                                    label_color = label_color,
                                    label_size = 10,
                                    label_angle = label_angle,
                                    label_position = "middle")


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


def get_gbk_blast_act(candidates_file, gbk_file, flanking_length, name_to_group_number_dict, path_to_output_act_folder):
    matches = open(candidates_file)
    total = 0
    for match in matches:
        total += 1

    matches = open(candidates_file)
    n = 1
    for match in matches:
        match = match.strip()
        #genes = sorted(match.split('\t'))
        genes = match.split('\t')
        folder_name = '___'.join(genes)
        stdout.write("\rProcessing %dth of %d HGT candidates: %s" % (n, total, folder_name))
        os.mkdir('%s/%s' % (path_to_output_act_folder, folder_name))

        # Extract gbk and faa files
        records = SeqIO.parse(gbk_file, 'genbank')
        gene_gc_dict = {}
        for record in records:
            for gene_f in record.features :
                if 'locus_tag' in gene_f.qualifiers:
                    for gene_1 in genes :
                        if gene_1 in gene_f.qualifiers["locus_tag"]:
                            start = gene_f.location.start.position
                            end = gene_f.location.end.position
                            orf = record.seq[start :end]
                            gene_gc = GC(orf)
                            gene_gc_dict[gene_1] = float("{0:.2f}".format(gene_gc))
                            pwd_gbk_file = '%s/%s/%s.gbk' % (path_to_output_act_folder, folder_name, gene_1)
                            pwd_fasta_file = '%s/%s/%s.fasta' % (path_to_output_act_folder, folder_name, gene_1)
                            SeqIO.write(record, pwd_gbk_file, 'genbank')
                            SeqIO.write(record, pwd_fasta_file, 'fasta')

                            # get flanking regions
                            get_flanking_region(pwd_gbk_file, gene_1, flanking_length)

        # Run Blast
        prefix_c = '%s/%s' % (path_to_output_act_folder, folder_name)
        output_c = '%s/%s.txt' % (prefix_c, folder_name)
        query_c = '%s/%s_%sbp.fasta' % (prefix_c, genes[0], flanking_length)
        subject_c = '%s/%s_%sbp.fasta' % (prefix_c, genes[1], flanking_length)

        parameters_c_n = '-evalue 1e-5 -outfmt 6 -task blastn'
        command_blast = 'blastn -query %s -subject %s -out %s %s' % (query_c, subject_c, output_c, parameters_c_n)
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
                                                            name = gene1_contig.name,
                                                            greytrack = True,
                                                            greytrack_labels = 1,
                                                            greytrack_font = 'Helvetica',
                                                            greytrack_fontsize = 12,
                                                            height = 0.35,
                                                            start = 0,
                                                            end = len(gene1_contig),
                                                            scale = True,
                                                            scale_fontsize = 6,
                                                            scale_ticks = 1,
                                                            scale_smalltick_interval = 10000,
                                                            scale_largetick_interval = 10000)
            # add gene content track for gene2_contig
            contig_2_gene_content_track = diagram.new_track(1,
                                                            name = gene2_contig.name,
                                                            greytrack = True,
                                                            greytrack_labels = 1,
                                                            greytrack_font = 'Helvetica',
                                                            greytrack_fontsize = 12,
                                                            height = 0.35,
                                                            start = 0,
                                                            end = len(gene2_contig),
                                                            scale = True,
                                                            scale_fontsize = 6,
                                                            scale_ticks = 1,
                                                            scale_smalltick_interval = 10000,
                                                            scale_largetick_interval = 10000)

            # add blank feature/graph sets to each track
            feature_sets_1 = contig_1_gene_content_track.new_set(type = 'feature')
            feature_sets_2 = contig_2_gene_content_track.new_set(type = 'feature')

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
                target = each_line_split[1]
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

            diagram.draw(format = "linear",
                         orientation = "landscape",
                         pagesize = (75 * cm, 25 * cm),
                         fragments = 1,
                         start = 0,
                         end = max_len)

            diagram.write('%s/%s.eps' % (path_to_output_act_folder, folder_name), "EPS")
            # shutil.rmtree('%s/%s' % (path_to_output_act_folder, folder_name), ignore_errors=True)
            # if os.path.isdir('%s/%s' % (path_to_output_act_folder, folder_name)):
            #     shutil.rmtree('%s/%s' % (path_to_output_act_folder, folder_name), ignore_errors=True)
        n += 1


def add_direction(input_file, output_file):
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
        each_split_0 = each.split('\t')[0]
        each_split_1 = each.split('\t')[1]
        each_split_0_bin = each_split_0.split('_')[0]
        each_split_1_bin = each_split_1.split('_')[0]
        output.write('%s\t%s<-%s\n' % (each, each_split_0_bin, each_split_1_bin))

    for each in overlap_list:
        output.write('%s\tN/A\n' % each)

    output.close()


############################################## Read in configuration file ##############################################

parser = argparse.ArgumentParser()
config = configparser.ConfigParser()

parser.add_argument('-cfg',
                    required=True,
                    help='path to configuration file')

args = vars(parser.parse_args())
pwd_cfg_file = args['cfg']

config.read(pwd_cfg_file)
grouping_file = config['FILES_AND_PARAMETERS']['grouping_file']
cover_cutoff = int(config['FILES_AND_PARAMETERS']['cover_cutoff'])
flanking_length = int(config['FILES_AND_PARAMETERS']['flanking_length'])
identity_percentile = int(config['FILES_AND_PARAMETERS']['identity_percentile'])
align_len_cutoff = int(config['FILES_AND_PARAMETERS']['align_len_cutoff'])

blast_results = config['FILES_AND_PARAMETERS']['blast_results']
gbk_file_name = config['FILES_AND_PARAMETERS']['gbk_file_name']

############################################### Define folder/file name ################################################

print('Define folder/file names and create output folder')

op_prefix = 'output_ip' + str(identity_percentile) + '%_al' + str(align_len_cutoff) + 'bp_c' + str(cover_cutoff) + '%'
op_prefix_iden_0 = 'output_al' + str(align_len_cutoff) + 'bp_c' + str(cover_cutoff) + '%'
op_folder = op_prefix
wd = os.getcwd()

iden_distrib_plot_folder =                          'Identity_distribution_images'
qual_idens_file =                                   'qualified_identities.txt'
qual_idens_file_gg =                                'qualified_identities_gg.txt'
qual_idens_file_gg_sorted =                         'qualified_identities_gg_sorted.txt'
unploted_groups_file =                              'unploted_groups.txt'
qual_idens_with_group_filename =                    'qualified_identities_with_group.txt'
qual_idens_with_group_sorted_filename =             'qualified_identities_with_group_sorted.txt'
subjects_in_one_line_filename =                     'subjects_in_one_line.txt'
group_pair_iden_cutoff_file_name =                  'Identity_cutoff.txt'
op_candidates_with_group_file_name =                'HGT_candidates_with_group.txt'
op_candidates_only_gene_file_name =                 'HGT_candidates.txt'
op_candidates_only_gene_file_name_with_direction =  'HGT_candidates_with_direction.txt'
gbk_subset_file =                                   'combined_subset.gbk'
op_act_folder_name =                                'ACT_images'

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
pwd_op_candidates_only_gene_file_with_direction =   '%s/%s/%s'    % (wd, op_folder, op_candidates_only_gene_file_name_with_direction)
pwd_gbk_subset_file =                               '%s/%s/%s'    % (wd, op_folder, gbk_subset_file)
pwd_op_act_folder =                                 '%s/%s/%s'    % (wd, op_folder, op_act_folder_name)
path_to_grouping_file =                             '%s/%s'       % (wd, grouping_file)
path_to_blast_results =                             '%s/%s'       % (wd, blast_results)
path_to_gbk_file =                                  '%s/%s'       % (wd, gbk_file_name)

########################################################################################################################

# forward to working directory
os.chdir(wd)

# create outputs folder
if os.path.isdir(op_folder):
    shutil.rmtree(op_folder, ignore_errors=True)
    if os.path.isdir(op_folder):
        shutil.rmtree(op_folder, ignore_errors=True)
    os.makedirs(op_folder)
else:
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

# plot overall identity distribution
all_identities_plot_title = 'All_vs_All'
plot_identity_list(all_identities, 'None', all_identities_plot_title, pwd_iden_distrib_plot_folder)

sleep(1.5)
print('\nPlotting identity distributions')

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
            # if choice in ['s', 'S']:
            #     pass
            # elif choice in ['p', 'P']:
            #     sleep(0.5)
            # restore current_group_pair_name and current_group_pair_identities for next group pair
            current_group_pair_name = group_pair
            current_group_pair_identities = []
            current_group_pair_identities.append(identity)
# for the last group
do()
group_pair_iden_cutoff_file.close()


sleep(1.5)
# if choice in ['p', 'P']:
#     print('\nAnalyzing Blast matches to get HGT candidates')
# else:
#     print('Analyzing Blast matches to get HGT candidates')
print('\nAnalyzing Blast matches to get HGT candidates')

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

# add direction to output file
add_direction(pwd_op_candidates_only_gene_file, pwd_op_candidates_only_gene_file_with_direction)


#################################################### Get ACT images ####################################################

sleep(1.5)
print('Preparing subset of combined gbk file for ACT plotting')

# get gene list of all candidates
all_candidates_genes = []
matches = open(pwd_op_candidates_only_gene_file)
for match in matches:
    match_split = match.strip().split('\t')
    for gene in match_split:
        if gene not in all_candidates_genes:
            all_candidates_genes.append(gene)


# get subset of combined gbk file
gbk_subset = open(pwd_gbk_subset_file, 'w')
records = SeqIO.parse(path_to_gbk_file, 'genbank')
records_recorded = []


for record in records:
    record_id = record.id
    for gene_f in record.features:
        if 'locus_tag' in gene_f.qualifiers:
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
gene_gc_content_dict = get_gbk_blast_act(pwd_op_candidates_only_gene_file, pwd_gbk_subset_file, flanking_length, name_to_group_number_dict, pwd_op_act_folder)


sleep(1.5)
print('\nRemove temporary files... ')
# remove temporary files
# os.remove(pwd_qual_iden_file)
# os.remove(pwd_qual_iden_file_gg)
# os.remove(pwd_qual_iden_file_gg_sorted)
# os.remove(pwd_qual_idens_with_group)
# os.remove(pwd_qual_idens_with_group_sorted)
# os.remove(pwd_subjects_in_one_line)
# os.remove(pwd_gbk_subset_file)

sleep(1.5)
print('All done, enjoy!')
