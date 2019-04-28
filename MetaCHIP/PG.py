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
import re
import glob
import shutil
import argparse
import warnings
import platform
import numpy as np
from ete3 import Tree
from Bio import SeqIO, AlignIO, Align
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
from MetaCHIP.MetaCHIP_config import config_dict


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


def get_group_num_from_grouping_file(grouping_file):
    group_list = set()
    for each_line in open(grouping_file):
        each_line_split = each_line.strip().split(',')
        group_list.add(each_line_split[0])

    group_num = len(group_list)

    return group_num


class BinRecord(object):

    def __init__(self, name, group, group_without_underscore):
        self.name = name
        self.group = group
        self.group_without_underscore = group_without_underscore


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def subset_tree(tree_file_in, leaf_node_list, tree_file_out):
    tree_in = Tree(tree_file_in, format=0)
    tree_in.prune(leaf_node_list, preserve_branch_length=True)
    tree_in.write(format=0, outfile=tree_file_out)


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


def get_number_of_group(grouping_file):

    group_list = []
    for each_genome in open(grouping_file):
        each_genome_split = each_genome.strip().split(',')
        group_id = each_genome_split[0]
        if group_id not in group_list:
            group_list.append(group_id)
    number_of_group = len(group_list)

    return number_of_group


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def get_species_tree_alignment(tmp_folder, path_to_prokka, path_to_hmm, pwd_hmmsearch_exe, pwd_mafft_exe):

    # Tests for presence of the tmp folder and deletes it
    if os.path.exists(tmp_folder):
        os.system('rm -r ' + tmp_folder)
    os.mkdir(tmp_folder)

    # List all prokka dirs in the target folder
    prokka_files = [i for i in os.listdir(path_to_prokka) if os.path.isdir(path_to_prokka + '/' + i)]
    print('Detected %i input genomes' % len(prokka_files))

    # Running hmmsearch on each file
    print('Running hmmsearch...')
    for f in prokka_files:
        # call hmmsearch
        os.system('%s -o /dev/null --domtblout %s/%s_hmmout.tbl %s %s/%s/%s.faa' % (pwd_hmmsearch_exe, tmp_folder, f, path_to_hmm, path_to_prokka, f, f))

        # Reading the protein file in a dictionary
        proteinSequence = {}
        for seq_record in SeqIO.parse('%s/%s/%s.faa' % (path_to_prokka, f, f), 'fasta'):
            proteinSequence[seq_record.id] = str(seq_record.seq)

        # Reading the hmmersearch table/extracting the protein part found beu hmmsearch out of the protein/Writing each protein sequence that was extracted to a fasta file (one for each hmm in phylo.hmm
        hmm_id = ''
        hmm_name = ''
        hmm_pos1 = 0
        hmm_pos2 = 0
        hmm_score = 0

        with open(tmp_folder + '/' + f.replace('prokka/', '') + '_hmmout.tbl', 'r') as tbl:
            for line in tbl:
                if line[0] == "#": continue
                line = re.sub('\s+', ' ', line)
                splitLine = line.split(' ')

                if (hmm_id == ''):
                    hmm_id = splitLine[4]
                    hmm_name = splitLine[0]
                    hmm_pos1 = int(splitLine[17]) - 1
                    hmm_pos2 = int(splitLine[18])
                    hmm_score = float(splitLine[13])
                elif (hmm_id == splitLine[4]):
                    if(float(splitLine[13]) > hmm_score):
                        hmm_name = splitLine[0]
                        hmm_pos1 = int(splitLine[17]) - 1
                        hmm_pos2 = int(splitLine[18])
                        hmm_score = float(splitLine[13])
                else:
                    file_out = open(tmp_folder + '/' + hmm_id + '.fasta', 'a+')
                    file_out.write('>' + f + '\n')
                    seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
                    file_out.write(str(seq) + '\n')
                    file_out.close()
                    hmm_id = splitLine[4]
                    hmm_name = splitLine[0]
                    hmm_pos1 = int(splitLine[17]) - 1
                    hmm_pos2 = int(splitLine[18])
                    hmm_score = float(splitLine[13])

            else:
                file_out = open(tmp_folder + '/' + hmm_id + '.fasta', 'a+')
                file_out.write('>' + f + '\n')
                seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
                file_out.write(str(seq) + '\n')
                file_out.close()

    # Call mafft to align all single fasta files with hmms
    files = os.listdir(tmp_folder)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    print('Running mafft...')
    for f in fastaFiles:
        fastaFile1 = '%s/%s' % (tmp_folder, f)
        fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
        os.system(pwd_mafft_exe + ' --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)


def get_species_tree_newick(tmp_folder, each_subset, pwd_fasttree_exe, tree_folder, tree_name):
    # concatenating the single alignments
    concatAlignment = {}
    for element in each_subset:
        concatAlignment[element] = ''

    # Reading all single alignment files and append them to the concatenated alignment
    files = os.listdir(tmp_folder)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    for f in fastaFiles:
        fastaFile = tmp_folder + '/' + f
        proteinSequence = {}
        alignmentLength = 0
        for seq_record_2 in SeqIO.parse(fastaFile, 'fasta'):
            proteinName = seq_record_2.id
            proteinSequence[proteinName] = str(seq_record_2.seq)
            alignmentLength = len(proteinSequence[proteinName])

        for element in each_subset:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    pwd_alignment_file = '%s/%s.aln' % (tree_folder, tree_name)
    pwd_newick_file = '%s/%s.newick' % (tree_folder, tree_name)

    file_out = open(pwd_alignment_file, 'w')
    # for element in prokka_files:
    for element in each_subset:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()

    # calling fasttree for tree calculation
    os.system('%s -quiet %s > %s' % (pwd_fasttree_exe, pwd_alignment_file, pwd_newick_file))

    # draw species tree
    #species_tree = Phylo.read('phylogenticTree.phy', 'newick')
    #Phylo.draw_ascii(species_tree)


def remove_low_cov_and_consensus_columns(alignment_file_in, minimal_cov, min_consensus, alignment_file_out):

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

    def remove_columns_from_msa(alignment_in, cols_to_remove):

        # get 0 based index of all wanted columns
        cols_to_remove_0_base = [(i - 1) for i in cols_to_remove]
        aln_cols_index_all = list(range(alignment_in.get_alignment_length()))
        aln_cols_index_wanted = []
        for i in aln_cols_index_all:
            if i not in cols_to_remove_0_base:
                aln_cols_index_wanted.append(i)

        # get wanted alignment segments
        wanted_segments = list_to_segments(aln_cols_index_wanted)

        # create an empty Alignment object
        alignment_new = Align.MultipleSeqAlignment([])
        for sequence in alignment_in:
            new_seq_object = Seq('')
            new_seq_record = SeqRecord(new_seq_object)
            new_seq_record.id = sequence.id
            new_seq_record.description = sequence.description
            alignment_new.append(new_seq_record)

        # add wanted columns to empty Alignment object
        for segment in wanted_segments:

            # for single column segment
            if segment[0] == segment[1]:
                segment_value = alignment_in[:, segment[0]]

                m = 0
                for each_seq in alignment_new:
                    each_seq.seq = Seq(str(each_seq.seq) + segment_value[m])
                    m += 1

            # for multiple columns segment
            else:
                segment_value = alignment_in[:, (segment[0]):(segment[1] + 1)]
                alignment_new += segment_value

        return alignment_new

    def remove_low_cov_columns(alignment_in, min_cov_cutoff):

        # get columns with low coverage
        sequence_number = len(alignment_in)
        total_col_num = alignment_in.get_alignment_length()
        low_cov_columns = []
        n = 0
        while n < total_col_num:
            current_column = alignment_in[:, n]
            dash_number = current_column.count('-')
            gap_percent = (dash_number / sequence_number) * 100

            if gap_percent > min_cov_cutoff:
                low_cov_columns.append(n + 1)

            n += 1

        # remove identified columns
        alignment_new = remove_columns_from_msa(alignment_in, low_cov_columns)

        return alignment_new

    def remove_low_consensus_columns(alignment_in, min_css_cutoff):

        # get columns with low coverage
        sequence_number = len(alignment_in)
        total_col_num = alignment_in.get_alignment_length()
        low_css_columns = []
        n = 0
        while n < total_col_num:
            current_column = alignment_in[:, n]

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
            if most_abundant_aa_percent < min_css_cutoff:
                low_css_columns.append(n + 1)

            n += 1

        # remove identified columns
        alignment_new = remove_columns_from_msa(alignment_in, low_css_columns)

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


def get_ctg_match_cate_and_identity_distribution_plot(pwd_candidates_file_ET, pwd_plot_ctg_match_cate, pwd_iden_distribution_plot_BM, pwd_iden_distribution_plot_PG):

    # read in prediction results
    HGT_num_BM_normal = 0
    HGT_num_BM_at_end = 0
    HGT_num_BM_ctg_aln = 0
    HGT_num_PG_normal = 0
    HGT_num_PG_at_end = 0
    HGT_num_PG_ctg_aln = 0
    identity_list_BM_normal = []
    identity_list_BM_end_match = []
    identity_list_BM_full_length_match = []
    identity_list_PG_normal = []
    identity_list_PG_end_match = []
    identity_list_PG_full_length_match = []
    for each_HGT in open(pwd_candidates_file_ET):
        if not each_HGT.startswith('Gene_1'):
            each_HGT_split = each_HGT.strip().split('\t')
            identity = float(each_HGT_split[4])
            end_match = each_HGT_split[5]
            full_length_match = each_HGT_split[6]
            PG_validation = each_HGT_split[7]

            # get number for normal
            if (end_match == 'no') and (full_length_match == 'no'):
                HGT_num_BM_normal += 1
                identity_list_BM_normal.append(identity)
                if PG_validation != 'NA':
                    HGT_num_PG_normal += 1
                    identity_list_PG_normal.append(identity)

            # get number for at_end
            if (end_match == 'yes') and (full_length_match == 'no'):
                HGT_num_BM_at_end += 1
                identity_list_BM_end_match.append(identity)
                if PG_validation != 'NA':
                    HGT_num_PG_at_end += 1
                    identity_list_PG_end_match.append(identity)

            # get number for Ctg_align
            if (end_match == 'no') and (full_length_match == 'yes'):
                HGT_num_BM_ctg_aln += 1
                identity_list_BM_full_length_match.append(identity)
                if PG_validation != 'NA':
                    HGT_num_PG_ctg_aln += 1
                    identity_list_PG_full_length_match.append(identity)

    ################################################## plot at_ends_stat ###################################################

    n_groups = 2
    normal_list = (HGT_num_BM_normal, HGT_num_PG_normal)
    at_end_list = (HGT_num_BM_at_end, HGT_num_PG_at_end)
    ctg_match_list = (HGT_num_BM_ctg_aln, HGT_num_PG_ctg_aln)

    # create plot
    index = np.arange(n_groups)
    bar_width = 0.15

    plt.bar(index, normal_list, bar_width, alpha=0.4, color='g', label='Normal', align='center')
    plt.bar(index + bar_width, at_end_list, bar_width, alpha=0.4, color='orange', label='End Match', align='center')
    plt.bar(index + bar_width * 2, ctg_match_list, bar_width, alpha=0.4, color='r', label='Full Length Match', align='center')

    plt.ylabel('Number of predicted HGTs')
    plt.title('Location of predicted HGTs')
    plt.xticks(index + 1.5 * bar_width, ('Best-match', 'Phylogenetic'))
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    plt.tight_layout()
    plt.savefig(pwd_plot_ctg_match_cate, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()

    #################################### plot identity distribution of identified HGTs #####################################

    ########## for BM approach ##########
    num_bins = 50
    combined_list_BM = (identity_list_BM_normal, identity_list_BM_end_match, identity_list_BM_full_length_match)
    color_list = ['g', 'orange', 'r']
    label_list = ['Normal', 'End Match', 'Full Length Match']
    plt.hist(combined_list_BM, num_bins, alpha=0.6, normed=0, stacked=1, linewidth=0, color=color_list, label=label_list, rwidth=0.85)
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    total_HGT_num_BM = len(identity_list_BM_normal) + len(identity_list_BM_end_match) + len(identity_list_BM_full_length_match)
    plt.title('Identity distribution of identified %s HGTs' % total_HGT_num_BM)
    plt.xlabel('Identity (%)')
    plt.ylabel('Number of identified HGT')

    # Get plot
    plt.tight_layout()
    plt.savefig(pwd_iden_distribution_plot_BM, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


    ########## for PG approach ##########
    num_bins = 50
    combined_list_PG = (identity_list_PG_normal, identity_list_PG_end_match, identity_list_PG_full_length_match)
    color_list = ['g', 'orange', 'r']
    label_list = ['Normal', 'End Match', 'Full Length Match']
    plt.hist(combined_list_PG, num_bins, alpha=0.6, normed=0, stacked=1, linewidth=0, color=color_list, label=label_list, rwidth=0.85)
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    total_HGT_num_PG = len(identity_list_PG_normal) + len(identity_list_PG_end_match) + len(identity_list_PG_full_length_match)
    plt.title('Identity distribution of identified %s HGTs' % total_HGT_num_PG)
    plt.xlabel('Identity (%)')
    plt.ylabel('Number of identified HGT')

    # Get plot
    plt.tight_layout()
    plt.savefig(pwd_iden_distribution_plot_PG, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


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
    pwd_species_tree_newick = '%s/%s_species_tree.newick'             % (pwd_tree_folder, each_to_process_concate)

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

        # remove columns in alignment
        remove_low_cov_and_consensus_columns(pwd_seq_file_1st_aln, 50, 50, pwd_seq_file_2nd_aln)

        # run fasttree
        cmd_fasttree = '%s -quiet -wag %s > %s' % (pwd_fasttree_exe, pwd_seq_file_2nd_aln, pwd_gene_tree_newick)
        os.system(cmd_fasttree)

        # Get species tree
        subset_tree(pwd_SCG_tree_all, genome_subset, pwd_species_tree_newick)

        # remove temp files
        os.remove(self_seq)
        os.remove(gene_tree_seq)
        os.remove(pwd_seq_file_1st_aln)
        os.remove(pwd_seq_file_2nd_aln)
        if non_self_seq_num > 0:
            os.remove(non_self_seq)
            os.remove(blast_output)
            os.remove(blast_output_sorted)
            os.remove(gene_tree_seq_uniq)


def Ranger_worker(argument_list):
        each_paired_tree = argument_list[0]
        pwd_ranger_inputs_folder = argument_list[1]
        pwd_tree_folder = argument_list[2]
        pwd_ranger_exe = argument_list[3]
        pwd_ranger_outputs_folder = argument_list[4]

        # define Ranger-DTL input file name
        each_paired_tree_concate = '___'.join(each_paired_tree)
        pwd_species_tree_newick = '%s/%s_species_tree.newick' % (pwd_tree_folder, each_paired_tree_concate)
        pwd_gene_tree_newick =    '%s/%s_gene_tree.newick'    % (pwd_tree_folder, each_paired_tree_concate)
        #each_paired_tree_concate_short = '%s___%s' % (each_paired_tree[0].split('_')[-1], each_paired_tree[1].split('_')[-1])

        if (os.path.isfile(pwd_species_tree_newick) == True) and (os.path.isfile(pwd_gene_tree_newick) == True):

            ranger_inputs_file_name = each_paired_tree_concate + '.txt'
            ranger_outputs_file_name = each_paired_tree_concate + '_ranger_output.txt'
            # ranger_outputs_file_name_bootstrap = each_paired_tree_concate + '_ranger_bootstrap.txt'

            pwd_ranger_inputs = '%s/%s' % (pwd_ranger_inputs_folder, ranger_inputs_file_name)
            pwd_ranger_outputs = '%s/%s' % (pwd_ranger_outputs_folder, ranger_outputs_file_name)
            # pwd_ranger_outputs_bootstrap = '%s/%s' % (pwd_ranger_outputs_folder, ranger_outputs_file_name_bootstrap)

            # pwd_current_ranger_outputs_folder = '%s/%s' % (pwd_ranger_outputs_folder, each_paired_tree_concate_short)


            # read in species tree
            species_tree = Tree(pwd_species_tree_newick, format=0)
            species_tree.resolve_polytomy(recursive=True)  # solving multifurcations
            species_tree.convert_to_ultrametric() # for dated mode

            # read in gene tree
            gene_tree = Tree(pwd_gene_tree_newick, format=0)
            gene_tree.resolve_polytomy(recursive=True)  # solving multifurcations


            ################################################################################################################

            # change species tree leaf name for Ranger-DTL2, replace "_" with "XXXXX", then, replace "." with "SSSSS"
            for each_st_leaf in species_tree:
                each_st_leaf_name = each_st_leaf.name

                # replace '-' with 'XXXXX'
                if '_' in each_st_leaf_name:
                    each_st_leaf_name_no_Underline = 'XXXXX'.join(each_st_leaf_name.split('_'))
                else:
                    each_st_leaf_name_no_Underline = each_st_leaf_name

                # replace '.' with 'SSSSS'
                if '.' in each_st_leaf_name_no_Underline:
                    each_st_leaf_name_no_Underline_no_dot = 'SSSSS'.join(each_st_leaf_name_no_Underline.split('.'))
                else:
                    each_st_leaf_name_no_Underline_no_dot= each_st_leaf_name_no_Underline

                # rename species tree leaf name
                each_st_leaf.name = each_st_leaf_name_no_Underline_no_dot


            # change gene tree leaf name for Ranger-DTL2, replace "_" with "XXXXX", then, replace "." with "SSSSS"
            for each_gt_leaf in gene_tree:
                each_gt_leaf_name = each_gt_leaf.name

                # replace '-' with 'XXXXX'
                if '_' in each_gt_leaf_name:
                    each_gt_leaf_name_no_Underline = 'XXXXX'.join(each_gt_leaf_name.split('_')[:-1])
                else:
                    each_gt_leaf_name_no_Underline = each_gt_leaf_name

                # replace '.' with 'SSSSS'
                if '.' in each_gt_leaf_name_no_Underline:
                    each_gt_leaf_name_no_Underline_no_dot = 'SSSSS'.join(each_gt_leaf_name_no_Underline.split('.'))
                else:
                    each_gt_leaf_name_no_Underline_no_dot = each_gt_leaf_name_no_Underline

                # rename gene tree leaf name
                each_gt_leaf.name = each_gt_leaf_name_no_Underline_no_dot


            ################################################################################################################

            # write species tree and gene tree to Ranger-DTL input file
            ranger_inputs_file = open(pwd_ranger_inputs, 'w')

            # dated mode
            ranger_inputs_file.write('%s\n%s\n' % (species_tree.write(format=5), gene_tree.write(format=5)))
            ranger_inputs_file.close()

            # create ranger_outputs_folder
            #force_create_folder(pwd_current_ranger_outputs_folder)

            # run Ranger-DTL
            ranger_parameters = '-q -D 2 -T 3 -L 1'
            ranger_cmd = '%s %s -i %s -o %s' % (pwd_ranger_exe, ranger_parameters, pwd_ranger_inputs, pwd_ranger_outputs)
            os.system(ranger_cmd)

        # # run ranger with 100 bootstrap
        # ranger_bootstrap = 1
        # while ranger_bootstrap <= 100:
        #     ranger_outputs_bootstrap = '%s/%s_bootstrap%s' % (pwd_current_ranger_outputs_folder, each_paired_tree_concate, ranger_bootstrap)
        #     ranger_cmd_bootstrap = '%s %s -i %s -o %s' % (pwd_ranger_exe, ranger_parameters, pwd_ranger_inputs, ranger_outputs_bootstrap)
        #     os.system(ranger_cmd_bootstrap)
        #     ranger_bootstrap += 1
        #
        # # AggregateRanger_cmd
        # current_wd = os.getcwd()
        # os.chdir(pwd_ranger_outputs_folder)
        #
        # AggregateRanger_cmd = '%s %s/%s_bootstrap > %s' % (pwd_AggregateRanger_exe, each_paired_tree_concate_short, each_paired_tree_concate, ranger_outputs_file_name_bootstrap)
        # os.system(AggregateRanger_cmd)
        # os.system('rm -r %s' % each_paired_tree_concate_short)
        # os.chdir(current_wd)


def PG(args, config_dict):

    output_prefix =             args['p']
    grouping_level =            args['r']
    grouping_file =             args['g']
    cover_cutoff =              args['cov']
    align_len_cutoff =          args['al']
    flanking_length_kbp =       args['flk']
    identity_percentile =       args['ip']
    end_match_identity_cutoff = args['ei']
    num_threads =               args['t']
    keep_quiet =                args['quiet']

    # read in config file
    pwd_ranger_exe = config_dict['ranger_linux']
    if platform.system() == 'Darwin':
        pwd_ranger_exe = config_dict['ranger_mac']

    pwd_mafft_exe =     config_dict['mafft']
    pwd_fasttree_exe =  config_dict['fasttree']
    pwd_blastp_exe =    config_dict['blastp']
    circos_HGT_R =      config_dict['circos_HGT_R']

    warnings.filterwarnings("ignore")


    #################################### find matched grouping file if not provided  ###################################

    if grouping_level is None:
        grouping_level = 'x'

    MetaCHIP_wd =   '%s_MetaCHIP_wd'              % output_prefix
    pwd_log_file =  '%s/%s_%s_PG_%s.log'  % (MetaCHIP_wd, output_prefix, grouping_level, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))


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


    ############################################### Define folder/file name ################################################

    MetaCHIP_op_folder = '%s_%s%s_HGTs_ip%s_al%sbp_c%s_ei%sbp_f%skbp' % (output_prefix, grouping_level, group_num, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)

    genome_size_file_name =                             '%s_all_genome_size.txt'                      % (output_prefix)
    tree_folder =                                       '%s_%s%s_PG_tree_folder'                      % (output_prefix, grouping_level, group_num)
    ranger_inputs_folder_name =                         '%s_%s%s_PG_Ranger_input'                     % (output_prefix, grouping_level, group_num)
    ranger_outputs_folder_name =                        '%s_%s%s_PG_Ranger_output'                    % (output_prefix, grouping_level, group_num)
    candidates_file_name =                              '%s_%s%s_HGTs_BM.txt'                         % (output_prefix, grouping_level, group_num)
    candidates_seq_file_name =                          '%s_%s%s_HGTs_BM_nc.fasta'                    % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET =                           '%s_%s%s_HGTs_PG.txt'                         % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET_validated =                 '%s_%s%s_HGTs_PG_validated.txt'               % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET_validated_STAT_png =        '%s_%s%s_HGTs_PG_validated_stats.png'         % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET_validated_STAT_group_txt =  '%s_%s%s_HGTs_PG_validated_group_stats.txt'   % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET_validated_STAT_genome_txt = '%s_%s%s_HGTs_PG_validated_genome_stats.txt'  % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET_validated_fasta_nc =        '%s_%s%s_HGTs_PG_nc.fasta'                    % (output_prefix, grouping_level, group_num)
    candidates_file_name_ET_validated_fasta_aa =        '%s_%s%s_HGTs_PG_aa.fasta'                    % (output_prefix, grouping_level, group_num)
    flanking_region_plot_folder_name =                  '%s_%s%s_Flanking_region_plots'               % (output_prefix, grouping_level, group_num)
    newick_tree_file =                                  '%s_%s%s_species_tree.newick'                 % (output_prefix, grouping_level, group_num)
    combined_faa_file =                                 '%s_%s%s_combined_faa.fasta'                  % (output_prefix, grouping_level, group_num)
    grouping_id_to_taxon_file_name =                    '%s_%s%s_group_to_taxon.txt'                  % (output_prefix, grouping_level, group_num)
    grouping_file_with_id_filename =                    '%s_%s%s_grouping_with_id.txt'                % (output_prefix, grouping_level, group_num)
    combined_faa_file_subset =                          '%s_%s%s_combined_subset.faa'                 % (output_prefix, grouping_level, group_num)
    plot_identity_distribution_BM =                     '%s_%s%s_plot_HGT_identity_BM.png'            % (output_prefix, grouping_level, group_num)
    plot_identity_distribution_PG =                     '%s_%s%s_plot_HGT_identity_PG.png'            % (output_prefix, grouping_level, group_num)
    plot_at_ends_number =                               '%s_%s%s_plot_ctg_match_category.png'         % (output_prefix, grouping_level, group_num)
    plot_circos =                                       '%s_%s%s_plot_circos_PG.png'                  % (output_prefix, grouping_level, group_num)
    HGT_query_to_subjects_filename =                    '%s_%s%s_HGT_query_to_subjects.txt'           % (output_prefix, grouping_level, group_num)

    normal_folder_name =                                '1_Plots_normal'
    normal_folder_name_PG_validated =                   '1_Plots_normal_PG_validated'
    at_ends_folder_name =                               '2_Plots_end_match'
    at_ends_folder_name_PG_validated =                  '2_Plots_end_match_PG_validated'
    full_contig_match_folder_name =                     '3_Plots_full_length_match'
    full_contig_match_folder_name_PG_validated =        '3_Plots_full_length_match_PG_validated'

    pwd_MetaCHIP_op_folder =                            '%s/%s'                                       % (MetaCHIP_wd, MetaCHIP_op_folder)
    pwd_candidates_file =                               '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name)
    pwd_candidates_seq_file =                           '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_seq_file_name)
    pwd_candidates_file_ET =                            '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET)
    pwd_candidates_file_ET_validated =                  '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET_validated)
    pwd_candidates_file_ET_validated_STAT_png =         '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET_validated_STAT_png)
    pwd_candidates_file_ET_validated_STAT_group_txt =   '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET_validated_STAT_group_txt)
    pwd_candidates_file_ET_validated_STAT_genome_txt =  '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET_validated_STAT_genome_txt)
    pwd_candidates_file_ET_validated_fasta_nc =         '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET_validated_fasta_nc)
    pwd_candidates_file_ET_validated_fasta_aa =         '%s/%s'                                       % (pwd_MetaCHIP_op_folder, candidates_file_name_ET_validated_fasta_aa)
    pwd_plot_identity_distribution_BM =                 '%s/%s'                                       % (pwd_MetaCHIP_op_folder, plot_identity_distribution_BM)
    pwd_plot_identity_distribution_PG =                 '%s/%s'                                       % (pwd_MetaCHIP_op_folder, plot_identity_distribution_PG)
    pwd_plot_at_ends_number =                           '%s/%s'                                       % (pwd_MetaCHIP_op_folder, plot_at_ends_number)
    pwd_plot_circos =                                   '%s/%s'                                       % (pwd_MetaCHIP_op_folder, plot_circos)
    pwd_flanking_region_plot_folder =                   '%s/%s'                                       % (pwd_MetaCHIP_op_folder, flanking_region_plot_folder_name)
    pwd_1_normal_folder =                               '%s/%s'                                       % (pwd_flanking_region_plot_folder, normal_folder_name)
    pwd_1_normal_folder_PG_validated =                  '%s/%s'                                       % (pwd_flanking_region_plot_folder, normal_folder_name_PG_validated)
    pwd_2_at_ends_folder =                              '%s/%s'                                       % (pwd_flanking_region_plot_folder, at_ends_folder_name)
    pwd_2_at_ends_folder_PG_validated =                 '%s/%s'                                       % (pwd_flanking_region_plot_folder, at_ends_folder_name_PG_validated)
    pwd_3_full_contig_match_folder =                    '%s/%s'                                       % (pwd_flanking_region_plot_folder, full_contig_match_folder_name)
    pwd_3_full_contig_match_folder_PG_validated =       '%s/%s'                                       % (pwd_flanking_region_plot_folder, full_contig_match_folder_name_PG_validated)
    pwd_ranger_inputs_folder =                          '%s/%s'                                       % (pwd_MetaCHIP_op_folder, ranger_inputs_folder_name)
    pwd_ranger_outputs_folder =                         '%s/%s'                                       % (pwd_MetaCHIP_op_folder, ranger_outputs_folder_name)
    pwd_tree_folder =                                   '%s/%s'                                       % (pwd_MetaCHIP_op_folder, tree_folder)
    pwd_combined_faa_file =                             '%s/%s'                                       % (MetaCHIP_wd, combined_faa_file)
    pwd_combined_faa_file_subset =                      '%s/%s'                                       % (pwd_MetaCHIP_op_folder, combined_faa_file_subset)
    pwd_genome_size_file =                              '%s/%s'                                       % (MetaCHIP_wd, genome_size_file_name)
    pwd_newick_tree_file =                              '%s/%s'                                       % (MetaCHIP_wd, newick_tree_file)
    pwd_grouping_id_to_taxon_file =                     '%s/%s'                                       % (MetaCHIP_wd, grouping_id_to_taxon_file_name)
    pwd_grouping_file_with_id =                         '%s/%s/%s'                                    % (MetaCHIP_wd, MetaCHIP_op_folder, grouping_file_with_id_filename)
    pwd_HGT_query_to_subjects_file =                    '%s/%s/%s'                                    % (MetaCHIP_wd, MetaCHIP_op_folder, HGT_query_to_subjects_filename)


    ###################################### store ortholog information into dictionary ######################################

    # create folders
    force_create_folder(pwd_tree_folder)

    # get list of match pair list
    candidates_list = []
    candidates_list_genes = set()
    for match_group in open(pwd_candidates_file):
        if not match_group.startswith('Gene_1'):
            match_group_split = match_group.strip().split('\t')[:2]
            candidates_list.append(match_group_split)
            candidates_list_genes.add(match_group_split[0])
            candidates_list_genes.add(match_group_split[1])

    if candidates_list == []:
        report_and_log(('No HGT detected by BM approach, program exited!'), pwd_log_file, keep_quiet)
        exit()

    # for report and log
    report_and_log(('Get gene/genome member in gene/species tree for each BM predicted HGT'), pwd_log_file, keep_quiet)


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
    report_and_log(('Prepare subset of %s for building gene tree' % combined_faa_file), pwd_log_file, keep_quiet)

    # uniq gene id list
    gene_id_uniq_set = set()
    for each_gene in gene_id_overall:
        gene_id_uniq_set.add(each_gene)

    # prepare combined_ffn file subset to speed up
    pwd_combined_faa_file_subset_handle = open(pwd_combined_faa_file_subset, 'w')
    for each_gene in SeqIO.parse(pwd_combined_faa_file, 'fasta'):
        if each_gene.id in gene_id_uniq_set:
            pwd_combined_faa_file_subset_handle.write('>%s\n' % each_gene.id)
            pwd_combined_faa_file_subset_handle.write('%s\n' % str(each_gene.seq))
    pwd_combined_faa_file_subset_handle.close()


    ################################## Extract gene sequences, run mafft and fasttree ##################################

    # for report and log
    report_and_log(('Get species/gene tree for %s BM approach identified HGTs with %s cores' % (len(candidates_list), num_threads)), pwd_log_file, keep_quiet)

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
    report_and_log(('Running Ranger-DTL2 with dated mode'), pwd_log_file, keep_quiet)

    # put multiple arguments in list
    list_for_multiple_arguments_Ranger = []
    for each_paired_tree in candidates_list:
        list_for_multiple_arguments_Ranger.append([each_paired_tree, pwd_ranger_inputs_folder, pwd_tree_folder, pwd_ranger_exe, pwd_ranger_outputs_folder])

    pool = mp.Pool(processes=num_threads)
    pool.map(Ranger_worker, list_for_multiple_arguments_Ranger)
    pool.close()
    pool.join()


    ########################################### parse Ranger-DTL prediction result #########################################

    # for report and log
    report_and_log(('Parsing Ranger prediction results'), pwd_log_file, keep_quiet)

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
                        donor_p = '_'.join(donor_p.split('XXXXX'))
                        donor_p = '.'.join(donor_p.split('SSSSS'))
                        recipient_p = recipient.split('-->')[1][1:]
                        recipient_p = '_'.join(recipient_p.split('XXXXX'))
                        recipient_p = '.'.join(recipient_p.split('SSSSS'))
                        predicted_transfer = donor_p + '-->' + recipient_p
                        predicted_transfers.append(predicted_transfer)


            candidate_2_predictions_dict[each_ranger_prediction_concate] = predicted_transfers

            # get two possible transfer situation
            candidate_split_gene = each_ranger_prediction_concate.split('___')
            candidate_split_gene_only_genome = []
            for each_candidate in candidate_split_gene:
                each_candidate_genome = '_'.join(each_candidate.split('_')[:-1])
                candidate_split_gene_only_genome.append(each_candidate_genome)

            possible_hgt_1 = '%s-->%s' % (candidate_split_gene_only_genome[0], candidate_split_gene_only_genome[1])
            possible_hgt_2 = '%s-->%s' % (candidate_split_gene_only_genome[1], candidate_split_gene_only_genome[0])
            possible_hgts = [possible_hgt_1, possible_hgt_2]
            candidate_2_possible_direction_dict[each_ranger_prediction_concate] = possible_hgts


    #################################################### combine results ###################################################

    # for report and log
    report_and_log(('Add Ranger-DTL predicted direction to HGT_candidates.txt'), pwd_log_file, keep_quiet)

    # add results to output file of best blast match approach
    combined_output_handle = open(pwd_candidates_file_ET, 'w')
    combined_output_validated_handle = open(pwd_candidates_file_ET_validated, 'w')
    combined_output_validated_header = 'Gene_1\tGene_2\tGene_1_group\tGene_2_group\tIdentity\tend_match\tfull_length_match\tDirection\n'
    combined_output_validated_handle.write(combined_output_validated_header)
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
                combined_output_validated_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_id, donor_genome_id, identity, end_break, Ctg_align, validated_prediction))
            combined_output_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_id, donor_genome_id, identity, end_break, Ctg_align, validated_prediction))
    combined_output_handle.close()
    combined_output_validated_handle.close()

    # export sequence of validated candidates
    combined_output_validated_fasta_nc_handle = open(pwd_candidates_file_ET_validated_fasta_nc, 'w')
    combined_output_validated_fasta_aa_handle = open(pwd_candidates_file_ET_validated_fasta_aa, 'w')
    for each_candidate in SeqIO.parse(pwd_candidates_seq_file, 'fasta'):
        if each_candidate.id in validated_candidate_list:
            # output nc sequences
            SeqIO.write(each_candidate, combined_output_validated_fasta_nc_handle, 'fasta')
            # output aa sequences
            each_candidate_aa = each_candidate
            each_candidate_aa.seq = each_candidate_aa.seq.translate()
            SeqIO.write(each_candidate_aa, combined_output_validated_fasta_aa_handle, 'fasta')
    combined_output_validated_fasta_nc_handle.close()
    combined_output_validated_fasta_aa_handle.close()

    # for report and log
    report_and_log(('Done for Phylogenetic approach!'), pwd_log_file, keep_quiet)


    ###################################### separate PG validated flanking region plots #####################################

    # create folders
    force_create_folder(pwd_1_normal_folder_PG_validated)
    force_create_folder(pwd_2_at_ends_folder_PG_validated)
    force_create_folder(pwd_3_full_contig_match_folder_PG_validated)


    for PG_HGT in open(pwd_candidates_file_ET):
        if PG_HGT != combined_output_validated_header:
            PG_HGT_split = PG_HGT.strip().split('\t')
            gene_1 = PG_HGT_split[0]
            gene_2 = PG_HGT_split[1]
            At_ends = PG_HGT_split[5]
            Ctg_align = PG_HGT_split[6]
            Direction = PG_HGT_split[7]
            possible_file_name_1 = '%s___%s.SVG' % (gene_1, gene_2)
            possible_file_name_2 = '%s___%s.SVG' % (gene_2, gene_1)

            # normal
            action = 'cp'
            if (At_ends == 'no') and (Ctg_align == 'no') and (Direction != 'NA'):
                pwd_possible_file_1 = '%s/%s' % (pwd_1_normal_folder, possible_file_name_1)
                pwd_possible_file_2 = '%s/%s' % (pwd_1_normal_folder, possible_file_name_2)
                if os.path.isfile(pwd_possible_file_1) is True:
                    os.system('%s %s %s/' % (action, pwd_possible_file_1, pwd_1_normal_folder_PG_validated))
                if os.path.isfile(pwd_possible_file_2) is True:
                    os.system('%s %s %s/' % (action, pwd_possible_file_2, pwd_1_normal_folder_PG_validated))

            # full length match
            if (At_ends == 'no') and (Ctg_align == 'yes') and (Direction != 'NA'):
                pwd_possible_file_1 = '%s/%s' % (pwd_3_full_contig_match_folder, possible_file_name_1)
                pwd_possible_file_2 = '%s/%s' % (pwd_3_full_contig_match_folder, possible_file_name_2)
                if os.path.isfile(pwd_possible_file_1) is True:
                    os.system('%s %s %s/' % (action, pwd_possible_file_1, pwd_3_full_contig_match_folder_PG_validated))
                if os.path.isfile(pwd_possible_file_2) is True:
                    os.system('%s %s %s/' % (action, pwd_possible_file_2, pwd_3_full_contig_match_folder_PG_validated))

            # end match
            if (At_ends == 'yes') and (Ctg_align == 'no') and (Direction != 'NA'):
                pwd_possible_file_1 = '%s/%s' % (pwd_2_at_ends_folder, possible_file_name_1)
                pwd_possible_file_2 = '%s/%s' % (pwd_2_at_ends_folder, possible_file_name_2)
                if os.path.isfile(pwd_possible_file_1) is True:
                    os.system('%s %s %s/' % (action, pwd_possible_file_1, pwd_2_at_ends_folder_PG_validated))
                if os.path.isfile(pwd_possible_file_2) is True:
                    os.system('%s %s %s/' % (action, pwd_possible_file_2, pwd_2_at_ends_folder_PG_validated))


    ########################################################################################################################
    ######################################################## Get_plot ######################################################
    ########################################################################################################################

    # for report and log
    report_and_log(('Plot stats of identified HGTs'), pwd_log_file, keep_quiet)

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

                if PG_validation != 'NA':
                    identity_list_not_end_PG.append(identity)
                    HGT_num_PG_not_at_end += 1

            else: # at_end == 'yes':
                HGT_num_BM_at_end += 1
                if PG_validation != 'NA':
                    HGT_num_PG_at_end += 1


    ################################ plot contig match category and HGT identity distribution ##############################

    get_ctg_match_cate_and_identity_distribution_plot(pwd_candidates_file_ET, pwd_plot_at_ends_number, pwd_plot_identity_distribution_BM, pwd_plot_identity_distribution_PG)


    ################################ plot number of HGT detected from each genome and group ################################

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
    group_id_with_taxon = []
    if grouping_level != 'x':

        for each_group_2_taxon in open(pwd_grouping_id_to_taxon_file):
            each_group_2_taxon_split = each_group_2_taxon.strip().split(',')
            group_2_taxon_dict[each_group_2_taxon_split[0]] = each_group_2_taxon_split[1]

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

    if grouping_level == 'x':
        HGT_PG_STAT_handle.write('Group\tSize(Mbp)\tHGT\tHGT/Mbp\n')
    else:
        HGT_PG_STAT_handle.write('Group\tSize(Mbp)\tHGT\tHGT/Mbp\tTaxon\n')

    for each_g in group_list_uniq:

        if grouping_level == 'x':
            for_out = '%s\t%s\t%s\t%s\n' % (each_g, float("{0:.3f}".format(group_to_length_dict[each_g])), group_list_uniq_count[n], group_list_uniq_count_normalized[n])
        else:
            for_out = '%s\t%s\t%s\t%s\t%s\n' % (each_g, float("{0:.3f}".format(group_to_length_dict[each_g])), group_list_uniq_count[n], group_list_uniq_count_normalized[n], group_2_taxon_dict[each_g])

        HGT_PG_STAT_handle.write(for_out)
        n += 1
    HGT_PG_STAT_handle.close()


    # get color list for group level subplot
    color_list_according_group_uniq = []
    for each in color_list_according_group:
        if color_list_according_group_uniq == []:
            color_list_according_group_uniq.append(each)
        elif each != color_list_according_group_uniq[-1]:
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
    plt.xticks(x_range_genome, genome_list_according_group_with_group, rotation=315, fontsize=xticks_fontsize_genome, horizontalalignment='left')
    plt.ylabel('Number of HGT')

    # subplot 2
    plt.subplot(222)
    if grouping_level == 'x':
        plt.bar(x_range_group, group_list_uniq_count, tick_label=group_list_uniq, align='center', alpha=0.5, linewidth=0, color=color_list_according_group_uniq)
        plt.xticks(x_range_group, group_list_uniq, rotation=315, fontsize=xticks_fontsize_group,horizontalalignment='left')
    else:
        plt.bar(x_range_group, group_list_uniq_count, tick_label=group_id_with_taxon, align='center', alpha=0.5, linewidth=0, color=color_list_according_group_uniq)
        plt.xticks(x_range_group, group_id_with_taxon, rotation=315, fontsize=xticks_fontsize_group,horizontalalignment='left')

    # subplot 3
    plt.subplot(223)
    plt.bar(x_range_genome, num_list_according_group_norm, alpha=0.5, linewidth=0, color=color_list_according_group)
    plt.xticks(x_range_genome, genome_list_according_group_with_group, rotation=315, fontsize=xticks_fontsize_genome, horizontalalignment='left')
    plt.xlabel('Genome')
    plt.ylabel('Number of HGT / Mbp sequences')

    # subplot 4
    plt.subplot(224)
    if grouping_level == 'x':
        plt.bar(x_range_group, group_list_uniq_count_normalized, tick_label=group_list_uniq, align='center', alpha=0.5, linewidth=0, color=color_list_according_group_uniq)
        plt.xticks(x_range_group, group_list_uniq, rotation=315, fontsize=xticks_fontsize_group,horizontalalignment='left')
    else:
        plt.bar(x_range_group, group_list_uniq_count_normalized, tick_label=group_id_with_taxon, align='center', alpha=0.5, linewidth=0, color=color_list_according_group_uniq)
        plt.xticks(x_range_group, group_id_with_taxon, rotation=315, fontsize=xticks_fontsize_group,horizontalalignment='left')
    plt.xlabel('Group')

    # plot layout
    plt.subplots_adjust(wspace=0.2, top=0.95)

    # save plot
    plt.tight_layout()
    plt.savefig(pwd_candidates_file_ET_validated_STAT_png, dpi=300)


    ################################################### Get_circlize_plot ##################################################

    # for predicted HGTs:
    # 1. not end match
    # 2. not full length match
    # 3. PG validated

    pwd_cir_plot_t1 =              '%s/%s_%s%s_cir_plot_t1.txt'              % (pwd_MetaCHIP_op_folder, output_prefix, grouping_level, group_num)
    pwd_cir_plot_t1_sorted =       '%s/%s_%s%s_cir_plot_t1_sorted.txt'       % (pwd_MetaCHIP_op_folder, output_prefix, grouping_level, group_num)
    pwd_cir_plot_t1_sorted_count = '%s/%s_%s%s_cir_plot_t1_sorted_count.txt' % (pwd_MetaCHIP_op_folder, output_prefix, grouping_level, group_num)
    pwd_cir_plot_matrix_filename = '%s/%s_%s%s_cir_plot_matrix.csv'          % (pwd_MetaCHIP_op_folder, output_prefix, grouping_level, group_num)


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
            Direction = each_split[7]
            if Genome_1 not in name2id_dict:
                name2id_dict[Genome_1] = Genome_1_ID
            if Genome_2 not in name2id_dict:
                name2id_dict[Genome_2] = Genome_2_ID
            transfers.append(Direction)

    tmp1 = open(pwd_cir_plot_t1, 'w')
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
    os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, pwd_cir_plot_matrix_filename, pwd_plot_circos))

    # for report and log
    report_and_log(('Gene flow plot exported to: %s' % plot_circos), pwd_log_file, keep_quiet)


    ################################################### remove tmp files ###################################################

    # for report and log
    report_and_log(('Deleting temporary files'), pwd_log_file, keep_quiet)

    # remove tmp files
    os.remove(pwd_cir_plot_t1)
    os.remove(pwd_cir_plot_t1_sorted)
    os.remove(pwd_cir_plot_t1_sorted_count)
    #os.remove(pwd_cir_plot_matrix_filename)
    #os.remove(pwd_combined_faa_file)
    os.remove(pwd_combined_faa_file_subset)

    os.system('rm -r %s' % pwd_ranger_inputs_folder)
    os.system('rm -r %s' % pwd_ranger_outputs_folder)
    os.system('rm -r %s' % pwd_tree_folder)


    # for report and log
    report_and_log(('All done!'), pwd_log_file, keep_quiet)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    # arguments for PG approach
    parser.add_argument('-p',             required=True,  help='output prefix')
    parser.add_argument('-r',             required=False, default=None, help='grouping rank')
    parser.add_argument('-g',             required=False, help='grouping file')
    parser.add_argument('-cov',           required=False, type=int, default=75, help='coverage cutoff, default: 75')
    parser.add_argument('-al',            required=False, type=int, default=200, help='alignment length cutoff, default: 200')
    parser.add_argument('-flk',           required=False, type=int, default=10, help='the length of flanking sequences to plot (Kbp), default: 10')
    parser.add_argument('-ip',            required=False, type=int, default=90, help='identity percentile, default: 90')
    parser.add_argument('-ei',            required=False, type=float, default=90, help='end match identity cutoff, default: 95')
    parser.add_argument('-t',             required=False, type=int, default=1, help='number of threads, default: 1')
    parser.add_argument('-force',         required=False, action="store_true", help='overwrite previous results')
    parser.add_argument('-quiet',         required=False, action="store_true", help='Do not report progress')

    args = vars(parser.parse_args())

    PG(args, config_dict)

