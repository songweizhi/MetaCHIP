from __future__ import division

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
import sys
import glob
import shutil
import argparse
import itertools
import subprocess
from time import sleep
from datetime import datetime
from string import ascii_uppercase
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqFeature as SF
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp


def get_genome_length(genome_file):
    genome_len = 0
    for each_seq in SeqIO.parse(genome_file, 'fasta'):
        genome_len += len(each_seq.seq)
    return genome_len


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


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def export_aa_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.protein)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def prodigal_parser(seq_file, sco_file, prefix, output_folder):

    bin_ffn_file =     '%s.ffn' % prefix
    bin_faa_file =     '%s.faa' % prefix
    bin_gbk_file =     '%s.gbk' % prefix
    pwd_bin_ffn_file = '%s/%s'  % (output_folder, bin_ffn_file)
    pwd_bin_faa_file = '%s/%s'  % (output_folder, bin_faa_file)
    pwd_bin_gbk_file = '%s/%s'  % (output_folder, bin_gbk_file)

    # get sequence id list
    id_to_sequence_dict = {}
    sequence_id_list = []
    for each_seq in SeqIO.parse(seq_file, 'fasta'):
        id_to_sequence_dict[each_seq.id] = str(each_seq.seq)
        sequence_id_list.append(each_seq.id)


    # get sequence to cds dict and sequence to transl_table dict
    current_seq_id = ''
    current_transl_table = ''
    current_seq_csd_list = []
    seq_to_cds_dict = {}
    seq_to_transl_table_dict = {}
    for each_cds in open(sco_file):
        if each_cds.startswith('# Sequence Data'):

            # add to dict
            if current_seq_id != '':
                seq_to_cds_dict[current_seq_id] = current_seq_csd_list
                seq_to_transl_table_dict[current_seq_id] = current_transl_table

            # reset value
            current_seq_id = each_cds.strip().split('=')[-1][1:-1].split(' ')[0]
            current_transl_table = ''
            current_seq_csd_list = []

        elif each_cds.startswith('# Model Data'):
            current_transl_table = each_cds.strip().split(';')[-2].split('=')[-1]

        else:
            current_seq_csd_list.append('_'.join(each_cds.strip().split('_')[1:]))

    seq_to_cds_dict[current_seq_id] = current_seq_csd_list
    seq_to_transl_table_dict[current_seq_id] = current_transl_table


    bin_gbk_file_handle = open(pwd_bin_gbk_file, 'w')
    bin_ffn_file_handle = open(pwd_bin_ffn_file, 'w')
    bin_faa_file_handle = open(pwd_bin_faa_file, 'w')
    gene_index = 1
    for seq_id in sequence_id_list:

        # create SeqRecord
        current_sequence = Seq(id_to_sequence_dict[seq_id])
        current_SeqRecord = SeqRecord(current_sequence, id=seq_id)
        current_SeqRecord.seq.alphabet = generic_dna
        transl_table = seq_to_transl_table_dict[seq_id]

        # add SeqFeature to SeqRecord
        for cds in seq_to_cds_dict[seq_id]:

            # define locus_tag id
            locus_tag_id = '%s_%s' % (prefix, "{:0>5}".format(gene_index))

            # define FeatureLocation
            cds_split = cds.split('_')
            cds_start = SF.ExactPosition(int(cds_split[0]))
            cds_end = SF.ExactPosition(int(cds_split[1]))
            cds_strand = cds_split[2]
            current_strand = None
            if cds_strand == '+':
                current_strand = 1
            if cds_strand == '-':
                current_strand = -1
            current_feature_location = FeatureLocation(cds_start, cds_end, strand=current_strand)

            # get nc sequence
            sequence_nc = ''
            if cds_strand == '+':
                sequence_nc = id_to_sequence_dict[seq_id][cds_start-1:cds_end]
            if cds_strand == '-':
                sequence_nc = str(Seq(id_to_sequence_dict[seq_id][cds_start-1:cds_end], generic_dna).reverse_complement())

            # translate to aa sequence
            sequence_aa = str(SeqRecord(Seq(sequence_nc)).seq.translate(table=transl_table))

            # remove * at the end
            sequence_aa = sequence_aa[:-1]

            # export nc and aa sequences
            export_dna_record(sequence_nc, locus_tag_id, '', bin_ffn_file_handle)
            export_aa_record(sequence_aa, locus_tag_id, '', bin_faa_file_handle)

            # Define feature type
            current_feature_type = 'CDS'

            # Define feature qualifiers
            current_qualifiers_dict = {}
            current_qualifiers_dict['locus_tag'] = locus_tag_id
            current_qualifiers_dict['transl_table'] = transl_table
            current_qualifiers_dict['translation'] = sequence_aa

            # Create a SeqFeature
            current_feature = SeqFeature(current_feature_location, type=current_feature_type, qualifiers=current_qualifiers_dict)

            # Append Feature to SeqRecord
            current_SeqRecord.features.append(current_feature)
            gene_index += 1

        # export to gbk file
        SeqIO.write(current_SeqRecord, bin_gbk_file_handle, 'genbank')

    bin_gbk_file_handle.close()
    bin_ffn_file_handle.close()
    bin_faa_file_handle.close()


def prodigal_worker(argument_list):

    input_genome = argument_list[0]
    input_genome_folder = argument_list[1]
    pwd_prodigal_exe = argument_list[2]
    nonmeta_mode = argument_list[3]
    pwd_prodigal_output_folder = argument_list[4]
    pwd_ffn_folder = argument_list[5]
    pwd_faa_folder = argument_list[6]
    pwd_gbk_folder = argument_list[7]

    # prepare command (according to Prokka)
    input_genome_basename, input_genome_ext = os.path.splitext(input_genome)
    pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    pwd_output_sco = '%s/%s.sco' % (pwd_prodigal_output_folder, input_genome_basename)

    prodigal_cmd_meta = '%s -f sco -q -c -m -g 11 -p meta -i %s -o %s' % (
    pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)
    prodigal_cmd_nonmeta = '%s -f sco -q -c -m -g 11 -i %s -o %s' % (
    pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)

    if nonmeta_mode is True:
        prodigal_cmd = prodigal_cmd_nonmeta
    else:
        prodigal_cmd = prodigal_cmd_meta

    os.system(prodigal_cmd)

    # prepare ffn, faa and gbk files from prodigal output
    prodigal_parser(pwd_input_genome, pwd_output_sco, input_genome_basename, pwd_prodigal_output_folder)

    # move file to separate folders
    os.system('mv %s/%s.ffn %s' % (pwd_prodigal_output_folder, input_genome_basename, pwd_ffn_folder))
    os.system('mv %s/%s.faa %s' % (pwd_prodigal_output_folder, input_genome_basename, pwd_faa_folder))
    os.system('mv %s/%s.gbk %s' % (pwd_prodigal_output_folder, input_genome_basename, pwd_gbk_folder))


def hmmsearch_worker(argument_list):
    faa_file_basename = argument_list[0]
    output_prefix = argument_list[1]
    MetaCHIP_wd = argument_list[2]
    pwd_SCG_tree_wd = argument_list[3]
    pwd_hmmsearch_exe = argument_list[4]
    path_to_hmm = argument_list[5]

    # run hmmsearch
    pwd_faa_file = '%s/%s_faa_files/%s.faa' % (MetaCHIP_wd, output_prefix, faa_file_basename)
    os.system('%s -o /dev/null --domtblout %s/%s_hmmout.tbl %s %s' % (
    pwd_hmmsearch_exe, pwd_SCG_tree_wd, faa_file_basename, path_to_hmm, pwd_faa_file))

    # Reading the protein file in a dictionary
    proteinSequence = {}
    for seq_record in SeqIO.parse(pwd_faa_file, 'fasta'):
        proteinSequence[seq_record.id] = str(seq_record.seq)

    # Reading the hmmersearch table/extracting the protein part found beu hmmsearch out of the protein/Writing
    # each protein sequence that was extracted to a fasta file (one for each hmm in phylo.hmm
    hmm_id = ''
    hmm_name = ''
    hmm_pos1 = 0
    hmm_pos2 = 0
    hmm_score = 0

    with open(pwd_SCG_tree_wd + '/' + faa_file_basename + '_hmmout.tbl', 'r') as tbl:
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
                if (float(splitLine[13]) > hmm_score):
                    hmm_name = splitLine[0]
                    hmm_pos1 = int(splitLine[17]) - 1
                    hmm_pos2 = int(splitLine[18])
                    hmm_score = float(splitLine[13])
            else:
                file_out = open(pwd_SCG_tree_wd + '/' + hmm_id + '.fasta', 'a+')
                file_out.write('>' + faa_file_basename + '\n')
                if hmm_name != '':
                    seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
                file_out.write(str(seq) + '\n')
                file_out.close()
                hmm_id = splitLine[4]
                hmm_name = splitLine[0]
                hmm_pos1 = int(splitLine[17]) - 1
                hmm_pos2 = int(splitLine[18])
                hmm_score = float(splitLine[13])

        else:
            file_out = open(pwd_SCG_tree_wd + '/' + hmm_id + '.fasta', 'a+')
            file_out.write('>' + faa_file_basename + '\n')
            if hmm_name != '':
                seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
            file_out.write(str(seq) + '\n')
            file_out.close()


def mafft_worker(argument_list):
    faa_file_basename = argument_list[0]
    pwd_SCG_tree_wd = argument_list[1]
    pwd_mafft_exe = argument_list[2]

    fastaFile1 = '%s/%s' % (pwd_SCG_tree_wd, faa_file_basename)
    fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')

    # High speed
    mafft_cmd = '%s --quiet --retree 1 %s > %s ; rm %s' % (pwd_mafft_exe, fastaFile1, fastaFile2, fastaFile1)

    # High accuracy
    #mafft_cmd = '%s --quiet --maxiterate 1000 --globalpair %s > %s ; rm %s' % (pwd_mafft_exe, fastaFile1, fastaFile2, fastaFile1)


    os.system(mafft_cmd)


##################################################### CONFIGURATION ####################################################


parser = argparse.ArgumentParser()

parser.add_argument('-i', required=True, help='input genome folder')

parser.add_argument('-t', required=True, help='taxonomic classification')

parser.add_argument('-l', required=False, default='c', help='grouping rank')

parser.add_argument('-x', required=False, default='fasta', help='file extension')

parser.add_argument('-p', required=True, help='output prefix')

parser.add_argument('-nonmeta', required=False, action="store_true", help='annotation as non-metagenome assembled genomes (non-MAG)')

parser.add_argument('-threads', required=False, type=int, default=1, help='number of threads')

parser.add_argument('-grouping_only', required=False, action="store_true", help='run grouping only, deactivate gene calling and phylogenetic tree building')

parser.add_argument('-quiet', required=False, action="store_true", help='Do not report progress')

args = vars(parser.parse_args())
input_genome_folder = args['i']
GTDB_output_file = args['t']
grouping_level = args['l']
file_extension = args['x']
output_prefix = args['p']
num_threads = args['threads']
grouping_only = args['grouping_only']
keep_quiet = args['quiet']
nonmeta_mode = args['nonmeta']


# get path to current script
pwd_Get_clusters_script = sys.argv[0]
Get_clusters_script_path, file_name = os.path.split(pwd_Get_clusters_script)
path_to_hmm = '%s/phylo.hmm' % Get_clusters_script_path
add_group_to_tree_R = '%s/add_group_to_tree.R' % Get_clusters_script_path


# read in config file
pwd_cfg_file = '%s/config.txt' % Get_clusters_script_path
program_path_dict = get_program_path_dict(pwd_cfg_file)

pwd_makeblastdb_exe = program_path_dict['makeblastdb']
pwd_blastn_exe = program_path_dict['blastn']
pwd_prodigal_exe = program_path_dict['prodigal']
pwd_hmmsearch_exe = program_path_dict['hmmsearch']
pwd_mafft_exe = program_path_dict['mafft']
pwd_fasttree_exe = program_path_dict['fasttree']


# check whether input genome exist
input_genome_re = '%s/*.%s' % (input_genome_folder, file_extension)
input_genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_re)]
if input_genome_file_list == []:
    print('No input genome detected, program exited!')
    exit()

input_genome_num = len(input_genome_file_list)

time_format = '[%Y-%m-%d %H:%M:%S] '


############################################### prepare files and folders ##############################################

MetaCHIP_wd = '%s_MetaCHIP_wd' % output_prefix
prodigal_output_folder = '%s_prodigal_output'    % output_prefix
ffn_folder =             '%s_ffn_files'          % output_prefix
faa_folder =             '%s_faa_files'          % output_prefix
gbk_folder =             '%s_gbk_files'          % output_prefix
combined_ffn_file =      '%s_combined.ffn'       % output_prefix
combined_faa_file =      '%s_combined.faa'       % output_prefix
blast_db_folder =        '%s_blastdb'            % output_prefix
blast_results_file =     '%s_all_vs_all_ffn.tab' % output_prefix
log_file_name =          '%s_PrepIn_%s.log'      % (output_prefix, datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

pwd_prodigal_output_folder = '%s/%s' % (MetaCHIP_wd, prodigal_output_folder)
pwd_ffn_folder =             '%s/%s' % (MetaCHIP_wd, ffn_folder)
pwd_faa_folder =             '%s/%s' % (MetaCHIP_wd, faa_folder)
pwd_gbk_folder =             '%s/%s' % (MetaCHIP_wd, gbk_folder)
pwd_log_file_name =          '%s/%s' % (MetaCHIP_wd, log_file_name)
pwd_combined_ffn_file =      '%s/%s' % (MetaCHIP_wd, combined_ffn_file)
pwd_combined_faa_file =      '%s/%s' % (MetaCHIP_wd, combined_faa_file)
pwd_blast_db_folder =        '%s/%s' % (MetaCHIP_wd, blast_db_folder)
pwd_blast_results =          '%s/%s' % (MetaCHIP_wd, blast_results_file)



############################################## read GTDB output into dict  #############################################

# read GTDB output into dict
taxon_assignment_dict = {}
for each_genome in open(GTDB_output_file):
    each_split = each_genome.strip().split('\t')
    bin_name = each_split[0]

    assignment_full = []
    if len(each_split) == 1:
        assignment_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    elif (len(each_split) > 1) and (';' in each_split[1]):
        assignment = each_split[1].split(';')
        if len(assignment) == 7:
            assignment_full = assignment
        if len(assignment) == 6:
            assignment_full = assignment + ['s__']
        if len(assignment) == 5:
            assignment_full = assignment + ['g__', 's__']
        if len(assignment) == 4:
            assignment_full = assignment + ['f__', 'g__', 's__']
        if len(assignment) == 3:
            assignment_full = assignment + ['o__', 'f__', 'g__', 's__']
        if len(assignment) == 2:
            assignment_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

    elif (len(each_split) > 1) and (';' not in each_split[1]):
        assignment_full = [each_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

    # get dict
    taxon_assignment_dict[bin_name] = assignment_full


# get all identified taxon at defined ranks
rank_to_position_dict = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
specified_rank_pos = rank_to_position_dict[grouping_level]
identified_taxon_list = []
for each_TaxonAssign in taxon_assignment_dict:
    specified_rank_id = taxon_assignment_dict[each_TaxonAssign][specified_rank_pos]
    if specified_rank_id not in identified_taxon_list:
        identified_taxon_list.append(specified_rank_id)


# get the id of genomes assigned to each taxon at specified level
taxon_2_genome_dict = {}
for each_taxon in identified_taxon_list:

    genome_list = []
    for genome in taxon_assignment_dict:
        if taxon_assignment_dict[genome][specified_rank_pos] == each_taxon:
            genome_list.append(genome)

    taxon_2_genome_dict[each_taxon] = genome_list


# get the number of ignored genome
ignored_genome_id = '%s__' % grouping_level
ignored_genome_num = 0
if ignored_genome_id in taxon_2_genome_dict:
    ignored_genome_num = len(taxon_2_genome_dict[ignored_genome_id])


rank_abbre_dict = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}


######################################################### report #######################################################

# report group number
group_num = 0
if ignored_genome_id in taxon_2_genome_dict:
    group_num = len(taxon_2_genome_dict) - 1
else:
    group_num = len(taxon_2_genome_dict)


################################################### get grouping file ##################################################

# get group id list
group_index_list = get_group_index_list()

# define grouping file name
genome_size_file_name =                 '%s_genome_size.txt' % (output_prefix)
grouping_file_name =                    '%s_grouping_%s%s.txt'       % (output_prefix, grouping_level, group_num)
grouping_plot_name =                    '%s_grouping_%s%s.png'       % (output_prefix, grouping_level, group_num)
grouping_id_to_taxon_tmp_file_name =    '%s_group_to_taxon_tmp.txt' % (output_prefix)
grouping_id_to_taxon_file_name =        '%s_group_to_taxon_%s.txt'  % (output_prefix, grouping_level)
excluded_genome_file_name =             '%s_excluded_genomes_%s.txt'   % (output_prefix, grouping_level)


pwd_genome_size_file =              '%s/%s' % (MetaCHIP_wd, genome_size_file_name)
pwd_grouping_file =                 '%s/%s' % (MetaCHIP_wd, grouping_file_name)
pwd_grouping_plot =                 '%s/%s' % (MetaCHIP_wd, grouping_plot_name)
pwd_grouping_id_to_taxon_tmp_file = '%s/%s' % (MetaCHIP_wd, grouping_id_to_taxon_tmp_file_name)
pwd_grouping_id_to_taxon_file =     '%s/%s' % (MetaCHIP_wd, grouping_id_to_taxon_file_name)
pwd_excluded_genome_file =          '%s/%s' % (MetaCHIP_wd, excluded_genome_file_name)

n = 0
for each_taxon in taxon_2_genome_dict:
    if each_taxon != ignored_genome_id:
        group_id = group_index_list[n]
        for genome in taxon_2_genome_dict[each_taxon]:
            for_write_1 = '%s,%s\n' % (group_id, genome)
            for_write_2 = '%s,%s\n' % (group_id, each_taxon)
        n += 1


################################################## plot grouping stats #################################################

# plot the number of genomes in each group
group_id_all = []
group_id_uniq = []
for each_group_assignment in open(pwd_grouping_file):
    group_id = each_group_assignment.strip().split(',')[0]
    if group_id not in group_id_uniq:
        group_id_uniq.append(group_id)
    group_id_all.append(group_id)

group_id_uniq_sorted = sorted(group_id_uniq)


# read group_2_taxon into dict
group_2_taxon_dict = {}
for each_group_2_taxon in open(pwd_grouping_id_to_taxon_file):
    each_group_2_taxon_split = each_group_2_taxon.strip().split(',')
    group_2_taxon_dict[each_group_2_taxon_split[0]] = each_group_2_taxon_split[1]


group_id_with_taxon = []
for each_group in group_id_uniq_sorted:
    each_group_new = '(%s) %s' % (each_group, group_2_taxon_dict[each_group])
    group_id_with_taxon.append(each_group_new)


group_id_uniq_count = []
for each_id in group_id_uniq_sorted:
    group_id_uniq_count.append(group_id_all.count(each_id))


################################################### export genome size #################################################

if grouping_only == False:

    # get input genome list
    input_genome_re = '%s/*.%s' % (input_genome_folder, file_extension)
    input_genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_re)]


################################################# if not grouping only #################################################

if grouping_only == False:

    ######################################## run prodigal with multiprocessing #########################################

    # get input genome list
    input_genome_re = '%s/*.%s' % (input_genome_folder, file_extension)
    input_genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_re)]

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Running Prodigal with %s cores.\n' % num_threads)
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Running Prodigal with %s cores.' % num_threads)

    list_for_multiple_arguments_Prodigal = []
    for input_genome in input_genome_file_list:
        list_for_multiple_arguments_Prodigal.append([input_genome, input_genome_folder, pwd_prodigal_exe, nonmeta_mode, pwd_prodigal_output_folder, pwd_ffn_folder, pwd_faa_folder, pwd_gbk_folder])


    ################################################# get species tree #################################################

    # define file name
    SCG_tree_wd =                   '%s_get_SCG_tree_wd'        % output_prefix
    combined_alignment_file =       '%s_species_tree.aln'       % output_prefix
    newick_tree_file =              '%s_species_tree.newick'    % output_prefix
    pwd_SCG_tree_wd =               '%s/%s'                     % (MetaCHIP_wd, SCG_tree_wd)
    pwd_combined_alignment_file =   '%s/%s'                     % (MetaCHIP_wd, combined_alignment_file)
    pwd_newick_tree_file =          '%s/%s'                     % (MetaCHIP_wd, newick_tree_file)

    os.mkdir(pwd_SCG_tree_wd)


    ######################################## run hmmsearch with multiprocessing ########################################

    faa_file_re = '%s/*.faa' % pwd_faa_folder
    faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
    faa_file_list = sorted(faa_file_list)

    faa_file_basename_list = []
    for faa_file in faa_file_list:
        faa_file_basename, faa_file_extension = os.path.splitext(faa_file)
        faa_file_basename_list.append(faa_file_basename)

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Running hmmsearch with %s cores.\n' % num_threads)
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Running hmmsearch with %s cores.' % num_threads)


    list_for_multiple_arguments_hmmsearch = []
    for faa_file_basename in faa_file_basename_list:
        list_for_multiple_arguments_hmmsearch.append([faa_file_basename, output_prefix, MetaCHIP_wd, pwd_SCG_tree_wd, pwd_hmmsearch_exe, path_to_hmm])


    pool = mp.Pool(processes=num_threads)
    pool.map(hmmsearch_worker, list_for_multiple_arguments_hmmsearch)
    pool.close()
    pool.join()


    ########################################## run mafft with multiprocessing ##########################################

    # Call mafft to align all single fasta files with hmms
    files = os.listdir(pwd_SCG_tree_wd)
    fastaFiles = [i for i in files if i.endswith('.fasta')]

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Running mafft with %s cores.\n' % num_threads)
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Running mafft with %s cores.' % num_threads)

    list_for_multiple_arguments_mafft = []
    for faa_file_basename in fastaFiles:
        list_for_multiple_arguments_mafft.append([faa_file_basename, pwd_SCG_tree_wd, pwd_mafft_exe])

    # run mafft
    pool = mp.Pool(processes=num_threads)
    pool.map(mafft_worker, list_for_multiple_arguments_mafft)
    pool.close()
    pool.join()


    ############################################# Concatenating alignments #############################################

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Concatenating alignments...\n')
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Concatenating alignments...')

    # concatenating the single alignments
    concatAlignment = {}
    for element in faa_file_basename_list:
        concatAlignment[element] = ''

    # Reading all single alignment files and append them to the concatenated alignment
    files = os.listdir(pwd_SCG_tree_wd)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    for faa_file_basename in fastaFiles:
        fastaFile = pwd_SCG_tree_wd + '/' + faa_file_basename
        proteinSequence = {}
        alignmentLength = 0
        for seq_record_2 in SeqIO.parse(fastaFile, 'fasta'):
            proteinName = seq_record_2.id
            proteinSequence[proteinName] = str(seq_record_2.seq)
            alignmentLength = len(proteinSequence[proteinName])

        for element in faa_file_basename_list:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    file_out = open(pwd_combined_alignment_file, 'w')
    for element in faa_file_basename_list:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()


    ################################################### run fasttree ###################################################

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Running fasttree...\n')
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Running fasttree...')

    # calling fasttree for tree calculation
    os.system('%s -quiet %s > %s' % (pwd_fasttree_exe, pwd_combined_alignment_file, pwd_newick_tree_file))

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Species tree was exported to %s\n' % newick_tree_file)
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Species tree was exported to %s' % newick_tree_file)


    ############################################### run all vs all blastn ##############################################

    # for report and log
    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'Running blastn with %s cores, be patient...\n' % num_threads)
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'Running blastn with %s cores, be patient...' % num_threads)

    # combine ffn and faa files
    os.system('cat %s/*.ffn > %s' % (pwd_ffn_folder, pwd_combined_ffn_file))
    os.system('cat %s/*.faa > %s' % (pwd_faa_folder, pwd_combined_faa_file))

    # make blast db and run all vs all blastn
    if os.path.isdir(pwd_blast_db_folder):
        shutil.rmtree(pwd_blast_db_folder)
        if os.path.isdir(pwd_blast_db_folder):
            shutil.rmtree(pwd_blast_db_folder)
    os.mkdir(pwd_blast_db_folder)

    makeblastdb_cmd = '%s -in %s/%s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_blast_db_folder, combined_ffn_file)
    blast_parameters = '-evalue 1e-5 -num_threads %s -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn' % num_threads
    all_vs_all_blast_cmd = '%s -query %s -db %s/%s -out %s %s' % (pwd_blastn_exe, pwd_combined_ffn_file, pwd_blast_db_folder, combined_ffn_file, pwd_blast_results, blast_parameters)

    os.system('cp %s %s' % (pwd_combined_ffn_file, pwd_blast_db_folder))
    os.system(makeblastdb_cmd)
    os.system(all_vs_all_blast_cmd)


    ############################################### for report and log file ##############################################

    with open(pwd_log_file_name, 'a') as log_handle:
        log_handle.write(datetime.now().strftime(time_format) + 'All protein orthologs within the input genomes need to be obtained for PG approach, you can get it with the following commands:\n')
        log_handle.write(datetime.now().strftime(time_format) + '$ get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -n 16 -d %s -X\n' % gbk_folder)
    if keep_quiet == 0:
        print(datetime.now().strftime(time_format) + 'All protein orthologs within the input genomes need to be obtained for PG approach, you can get it with the following commands:')
        print(datetime.now().strftime(time_format) + '$ get_homologues.pl -f 70 -t 3 -S 70 -E 1e-05 -C 70 -G -n 16 -d %s -X' % gbk_folder)


