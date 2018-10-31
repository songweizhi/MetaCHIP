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
import warnings
import argparse
import itertools
import subprocess
import numpy as np
from ete3 import Tree
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
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram


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


def get_ffn_faa_from_gbk(input_gbk_folder, pwd_ffn_folder, pwd_faa_folder):

    # get input gbk file list
    input_gbk_re = '%s/*.gbk' % input_gbk_folder
    input_gbk_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_gbk_re)]

    # get ffn and faa file from input gbk file
    for gbk_file in input_gbk_file_list:

        # prepare file name
        gbk_file_basename, gbk_file_extension = os.path.splitext(gbk_file)
        pwd_gbk_file = '%s/%s' % (input_gbk_folder, gbk_file)
        pwd_output_ffn_file = '%s/%s.ffn' % (pwd_ffn_folder, gbk_file_basename)
        pwd_output_faa_file = '%s/%s.faa' % (pwd_faa_folder, gbk_file_basename)

        output_ffn_handle = open(pwd_output_ffn_file, 'w')
        output_faa_handle = open(pwd_output_faa_file, 'w')
        for seq_record in SeqIO.parse(pwd_gbk_file, 'genbank'):
            for feature in seq_record.features:
                if feature.type == "CDS":
                    seq_record_sequence = str(seq_record.seq)

                    # get DNA sequence
                    seq_nc = ''
                    if feature.location.strand == 1:
                        seq_nc = seq_record_sequence[feature.location.start:feature.location.end]
                    if feature.location.strand == -1:
                        nc_seq_rc = seq_record_sequence[feature.location.start:feature.location.end]
                        seq_nc = str(Seq(nc_seq_rc, generic_dna).reverse_complement())

                    # get aa sequence
                    #seq_aa = feature.qualifiers['translation'][0]
                    feature_id = feature.qualifiers['locus_tag'][0]
                    feature_description = feature.qualifiers['product'][0]

                    # export to file
                    export_dna_record(seq_nc, feature_id, feature_description, output_ffn_handle)
                    #export_aa_record(seq_aa, feature_id, feature_description, output_faa_handle)

        output_ffn_handle.close()
        output_faa_handle.close()


def get_rank_assignment_dict(rank, taxon_assignment_lineage_file):

    rank_to_position_dict = {'d': 1, 'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6, 's': 7}
    rank_position = rank_to_position_dict[rank]

    assignment_dict = {}
    for each in open(taxon_assignment_lineage_file):
        each_split = each.strip().split('\t')
        bin_name = each_split[0]
        assignment = each_split[1].split(';')
        assignment_no_num = []
        for each_assign in assignment:
            assignment_no_num.append(each_assign.split('(')[0])

        new_assign = ''
        if len(assignment_no_num) <= rank_position:
            new_assign = assignment_no_num[-1]
        elif len(assignment_no_num) > rank_position:
            new_assign = assignment_no_num[rank_position-1]

        assignment_dict[bin_name] = new_assign

    return assignment_dict


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


def get_node_distance(tree, node_1, node_2):
    distance = tree.get_distance(node_1, node_2)
    return distance


def get_group(n):
    leaf_name = leaf_node_list[n]
    grouping_id = bin_to_grouping_dict[leaf_name]
    return '(%s)_(%s)' % (grouping_id, leaf_name)


def get_taxon(n):
    leaf_name = leaf_node_list[n]
    taxon_assign = bin_to_taxon_dict[leaf_name]
    grouping_id = bin_to_grouping_dict[leaf_name]
    return '(%s)_(%s)_(%s)' % (grouping_id, leaf_name, taxon_assign)
    #return taxon_assign


def get_distance_matrix(tree_file):

    # read in tree
    tree_in = Tree(tree_file, format=0)

    # get leaf node list
    leaf_node_list = []
    for leaf_node in tree_in:
        leaf_node_list.append(leaf_node.name)
    leaf_node_list = sorted(leaf_node_list)

    # get list of distance list
    all_distances_lol = []
    for each_node_1 in leaf_node_list:
        current_node_distance_list = []
        for each_node_2 in leaf_node_list:
            distance = 0
            if each_node_1 != each_node_2:
                distance = get_node_distance(tree_in, each_node_1, each_node_2)
                distance = float("{0:.5f}".format(distance))
            current_node_distance_list.append(str(distance))
        all_distances_lol.append(current_node_distance_list)

    return all_distances_lol


def plot_clustering_dendrogram(cluster, leaf_font_size, leaf_label_func, color_threshold, pwd_png_file):

    plt.figure(figsize=(9, 15))
    plt.xlabel('Distance')
    dendrogram(cluster, orientation='left', leaf_rotation=0, leaf_font_size=leaf_font_size, leaf_label_func=leaf_label_func, color_threshold=color_threshold)
    plt.axvline(x=max_d, c='k', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(pwd_png_file, dpi=300)
    plt.close()


##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()

parser.add_argument('-i', required=True, help='input genome folder')

parser.add_argument('-x', required=True, help='file extension')

parser.add_argument('-p', required=True, help='output prefix')

parser.add_argument('-dc', required=False, type=float, default=None, help='distance cutoff')

parser.add_argument('-fs', required=False, type=int, default=9, help='leaf name font size')

parser.add_argument('-ls', required=False, type=float, default=0.01, help='label shift on the tree plot')

parser.add_argument('-taxon', required=False, default=None, help='taxonomy classification of input genomes, if available')

parser.add_argument('-tr', required=False, default='c', help='taxon ranks')

args = vars(parser.parse_args())

input_genome_folder = args['i']
file_extension = args['x']
output_prefix = args['p']
max_d = args['dc']
leaf_font_size = args['fs']
taxon_classification_file = args['taxon']
taxon_rank = args['tr']
label_shift = args['ls']

# get path to current script
pwd_Get_clusters_script = sys.argv[0]
Get_clusters_script_path, file_name = os.path.split(pwd_Get_clusters_script)
path_to_hmm = '%s/phylo.hmm' % Get_clusters_script_path
add_group_to_tree_R = '%s/add_group_to_tree.R' % Get_clusters_script_path

# read in config file
pwd_cfg_file = '%s/config.txt' % Get_clusters_script_path
program_path_dict = get_program_path_dict(pwd_cfg_file)
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

########################################################################################################################

# disable warnings
warnings.filterwarnings("ignore")

# create MetaCHIP outputs folder
MetaCHIP_wd = '%s_MetaCHIP_wd' % output_prefix

prodigal_output_folder = '%s_prodigal_output' % output_prefix
ffn_folder =             '%s_ffn_files'       % output_prefix
faa_folder =             '%s_faa_files'       % output_prefix
gbk_folder =             '%s_gbk_files'       % output_prefix

pwd_prodigal_output_folder = '%s/%s' % (MetaCHIP_wd, prodigal_output_folder)
pwd_ffn_folder =             '%s/%s' % (MetaCHIP_wd, ffn_folder)
pwd_faa_folder =             '%s/%s' % (MetaCHIP_wd, faa_folder)
pwd_gbk_folder =             '%s/%s' % (MetaCHIP_wd, gbk_folder)

if os.path.isdir(MetaCHIP_wd):
    shutil.rmtree(MetaCHIP_wd, ignore_errors=True)
    if os.path.isdir(MetaCHIP_wd):
        shutil.rmtree(MetaCHIP_wd, ignore_errors=True)
        if os.path.isdir(MetaCHIP_wd):
            shutil.rmtree(MetaCHIP_wd, ignore_errors=True)
    os.mkdir(MetaCHIP_wd)
    os.mkdir(pwd_prodigal_output_folder)
    os.mkdir(pwd_ffn_folder)
    os.mkdir(pwd_faa_folder)
    os.mkdir(pwd_gbk_folder)
else:
    os.mkdir(MetaCHIP_wd)
    os.mkdir(pwd_prodigal_output_folder)
    os.mkdir(pwd_ffn_folder)
    os.mkdir(pwd_faa_folder)
    os.mkdir(pwd_gbk_folder)

##################################################### run prodigal #####################################################

# get input genome list
input_genome_re = '%s/*.%s' % (input_genome_folder, file_extension)
input_genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_re)]

# report current processing
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running Prodigal.')


MaxProcesses = 10
Processes = []
for input_genome in input_genome_file_list:

    # prepare command (according to Prokka)
    input_genome_basename, input_genome_ext = os.path.splitext(input_genome)
    pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    pwd_output_sco = '%s/%s.sco' % (pwd_prodigal_output_folder, input_genome_basename)
    prodigal_cmd = '%s -f sco -q -c -m -g 11 -p meta -i %s -o %s' % (pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)


    # run with subprocess
    prodigal_cmd_list = [pwd_prodigal_exe, '-f', 'sco',  '-q', '-c', '-m', '-g', '11', '-p', 'meta', '-i', pwd_input_genome, '-o', pwd_output_sco]

    # keep wait if there is no spare slots
    while len(Processes) >= MaxProcesses:
        sleep(0.1)
        for process in Processes:
            if process.poll() is not None:
                Processes.remove(process)

    # submit new subprocess
    p = subprocess.Popen(prodigal_cmd_list)
    Processes.append(p)

# wait for completion
for p in Processes:
    p.wait()


for input_genome_2 in input_genome_file_list:

    # prepare filename
    input_genome_basename, input_genome_ext = os.path.splitext(input_genome_2)
    pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome_2)
    pwd_output_sco = '%s/%s.sco' % (pwd_prodigal_output_folder, input_genome_basename)

    # prepare ffn, faa and gbk files from prodigal output
    prodigal_parser(pwd_input_genome, pwd_output_sco, input_genome_basename, pwd_prodigal_output_folder)

    # move file to separate folders
    os.system('mv %s/%s.ffn %s' % (pwd_prodigal_output_folder, input_genome_basename, pwd_ffn_folder))
    os.system('mv %s/%s.faa %s' % (pwd_prodigal_output_folder, input_genome_basename, pwd_faa_folder))
    os.system('mv %s/%s.gbk %s' % (pwd_prodigal_output_folder, input_genome_basename, pwd_gbk_folder))


################################################### get species tree ###################################################

# define file name
SCG_tree_wd =                   '%s_get_SCG_tree_wd'        % output_prefix
combined_alignment_file =       '%s_species_tree.aln'       % output_prefix
newick_tree_file =              '%s_species_tree.newick'    % output_prefix
pwd_SCG_tree_wd =               '%s/%s'                     % (MetaCHIP_wd, SCG_tree_wd)
pwd_combined_alignment_file =   '%s/%s'                     % (MetaCHIP_wd, combined_alignment_file)
pwd_newick_tree_file =          '%s/%s'                     % (MetaCHIP_wd, newick_tree_file)
os.mkdir(pwd_SCG_tree_wd)


faa_file_re = '%s/*.faa' % pwd_faa_folder
faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
faa_file_list = sorted(faa_file_list)

faa_file_basename_list = []
for faa_file in faa_file_list:
    faa_file_basename, faa_file_extension = os.path.splitext(faa_file)
    faa_file_basename_list.append(faa_file_basename)


# report current processing
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running hmmsearch.')

for faa_file_basename in faa_file_basename_list:

    # run hmmsearch
    pwd_faa_file = '%s/%s_faa_files/%s.faa' % (MetaCHIP_wd, output_prefix, faa_file_basename)
    os.system('%s -o /dev/null --domtblout %s/%s_hmmout.tbl %s %s' % (pwd_hmmsearch_exe, pwd_SCG_tree_wd, faa_file_basename, path_to_hmm, pwd_faa_file))

    # Reading the protein file in a dictionary
    proteinSequence = {}
    for seq_record in SeqIO.parse(pwd_faa_file, 'fasta'):
        proteinSequence[seq_record.id] = str(seq_record.seq)

    # Reading the hmmersearch table/extracting the protein part found beu hmmsearch out of the protein/Writing each protein sequence that was extracted to a fasta file (one for each hmm in phylo.hmm
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

# Call mafft to align all single fasta files with hmms
files = os.listdir(pwd_SCG_tree_wd)
fastaFiles = [i for i in files if i.endswith('.fasta')]
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running mafft...')
for faa_file_basename in fastaFiles:
    fastaFile1 = '%s/%s' % (pwd_SCG_tree_wd, faa_file_basename)
    fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
    os.system(pwd_mafft_exe + ' --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)


# concatenating the single alignments
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Concatenating alignments...')
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

# calling fasttree for tree calculation
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running fasttree...')
os.system('%s -quiet %s > %s' % (pwd_fasttree_exe, pwd_combined_alignment_file, pwd_newick_tree_file))

# Decomment the two following lines if tree is rooted but should be unrooted
# phyloTree = dendropy.Tree.get(path='phylogenticTree.phy', schema='newick', rooting='force-unrooted')
# dendropy.Tree.write_to_path(phyloTree, 'phylogenticTree_unrooted.phy', 'newick')
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' The built species tree was exported to %s' % newick_tree_file)

###################################################### get cluster #####################################################

# get bin to taxon dict
bin_to_taxon_dict = {}
if taxon_classification_file != None:
    bin_to_taxon_dict = get_rank_assignment_dict(taxon_rank, taxon_classification_file)

# read in tree
tree_in = Tree(pwd_newick_tree_file, format=0)

# get sorted leaf node list
leaf_node_list = []
for leaf_node in tree_in:
    leaf_node_list.append(leaf_node.name)
leaf_node_list = sorted(leaf_node_list)

# get distance matrix from tree file
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get distance matrix from tree')

all_distances_lol = get_distance_matrix(pwd_newick_tree_file)

# turn list of distance list into arrary
all_distances_lol_array = np.array(all_distances_lol)

# get linkage
cluster = linkage(all_distances_lol_array, method='single')

# get maximum distance for clustering
distance_list = []
for each in cluster:
    distance_list.append(each[2])

# get distance cutoff
percentile_for_distances_cutoff = 90
distance_of_percentile = np.percentile(distance_list, percentile_for_distances_cutoff)

if max_d == None:
    max_d = distance_of_percentile
    sleep(0.5)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Determined distance cutoff is: %s, you can change it with option "-dc"' % float("{0:.2f}".format(max_d)))
else:
    sleep(0.5)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Distance cutoff specified to %s' % max_d)

# get flat clusters
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping input genomes based on above distance cutoff')
flat_clusters = fcluster(cluster, max_d, criterion='distance')

# get group number
group_index_list = []
for each_index in flat_clusters:
    if int(each_index) not in group_index_list:
        group_index_list.append(each_index)
group_number = len(group_index_list)

# define output file name with grouping included
png_file_group =            '%s_grouping_g%s.png'       % (output_prefix, group_number)
grouping_file =             '%s_grouping_g%s.txt'       % (output_prefix, group_number)
grouping_file_temp =        '%s_grouping_g%s_tmp.file'  % (output_prefix, group_number)
pwd_png_file_group =        '%s/%s'                     % (MetaCHIP_wd, png_file_group)
pwd_grouping_file =         '%s/%s'                     % (MetaCHIP_wd, grouping_file)
pwd_grouping_file_temp =    '%s/%s'                     % (MetaCHIP_wd, grouping_file_temp)

# get grouping file
group_index_list = get_group_index_list()
grouping_file_temp_handle = open(pwd_grouping_file_temp, 'w')
bin_to_grouping_dict = {}
n = 0
for each_leaf in leaf_node_list:
    leaf_cluster_index = int(flat_clusters[n])
    leaf_grouping_id = group_index_list[leaf_cluster_index - 1]
    grouping_file_temp_handle.write('%s,%s\n' % (leaf_grouping_id,each_leaf))
    bin_to_grouping_dict[each_leaf] = leaf_grouping_id
    n += 1
grouping_file_temp_handle.close()

# sort grouping file
os.system('cat %s | sort > %s; rm %s' % (pwd_grouping_file_temp, pwd_grouping_file, pwd_grouping_file_temp))

# report
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping results exported to: %s' % grouping_file)

# calculate full dendrogram
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get plot for visualization')
if taxon_classification_file != None:
    plot_clustering_dendrogram(cluster, leaf_font_size, get_taxon, max_d, pwd_png_file_group)
else:
    plot_clustering_dendrogram(cluster, leaf_font_size, get_group, max_d, pwd_png_file_group)

# call R
current_wd = os.getcwd()
os.chdir(MetaCHIP_wd)
add_group_to_tree_R_cmd = 'Rscript %s -t %s -g %s -l %s > /dev/null' % (add_group_to_tree_R, newick_tree_file, grouping_file, label_shift)
print(add_group_to_tree_R_cmd)
os.system(add_group_to_tree_R_cmd)
os.chdir(current_wd)

# report done
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping step done!')

sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' You may want to modify the grouping results based on the taxonomy\nof your input genomes/bins, you can do this by changing their group assignment specified\nin the first column of %s' % grouping_file)

