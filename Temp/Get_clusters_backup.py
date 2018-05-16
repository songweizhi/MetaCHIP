#!/usr/bin/python
import re
import os
import warnings
import argparse
import itertools
import numpy as np
from ete3 import Tree
from time import sleep
from datetime import datetime
from string import ascii_uppercase
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram


# to do
# selected_rank = 'c'
# output_prefix = 'NorthSea'
# add_group_to_tree_R = '~/R_scripts/newick_tree/add_group_to_tree.R'


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def get_bin_to_taxon_dict(taxon_classification_file):
    bin_to_taxon_dict = {}
    for each_assignment in open(taxon_classification_file):
        each_assignment_split = each_assignment.strip().split('\t')
        bin_name = each_assignment_split[0]
        taxon_assign = each_assignment_split[1]
        bin_to_taxon_dict[bin_name] = taxon_assign
    return bin_to_taxon_dict


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
    tree_in = Tree(tree_file, format=3)

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


######################################################### input ########################################################

parser = argparse.ArgumentParser()

parser.add_argument('-prokka_output',
                    required=True,
                    help='the folder holds Prokka outputs for input genomes')

parser.add_argument('-hmm',
                    required=True,
                    help='the phylo.hmm file')

parser.add_argument('-dc',
                    required=False,
                    type=float,
                    default=None,
                    help='distance cutoff')

parser.add_argument('-fs',
                    required=False,
                    type=int,
                    default=6,
                    help='leaf name font size')

parser.add_argument('-taxon',
                    required=False,
                    default=None,
                    help='taxon classification if available')

parser.add_argument('-hmmsearch',
                    required=False,
                    default='hmmsearch',
                    help='path to hmmsearch executable')

parser.add_argument('-mafft',
                    required=False,
                    default='mafft',
                    help='path to mafft executable')

parser.add_argument('-fasttree',
                    required=False,
                    default='fasttree',
                    help='path to fasttree executable')

args = vars(parser.parse_args())
path_to_prokka = args['prokka_output']
path_to_hmm = args['hmm']
pwd_hmmsearch_exe = args['hmmsearch']
pwd_mafft_exe = args['mafft']
pwd_fasttree_exe = args['fasttree']
max_d = args['dc']
leaf_font_size = args['fs']
taxon_classification_file = args['taxon']
selected_rank = 'c'
output_prefix = 'NorthSea'
add_group_to_tree_R = '~/R_scripts/newick_tree/add_group_to_tree.R'

#################################################### define file name ##################################################

# disable warnings
warnings.filterwarnings("ignore")
wd = os.getcwd()

tmp_folder =                '%s_get_SCG_tree_tmp'       % output_prefix
combined_alignment_file =   '%s_species_tree.aln'       % output_prefix
newick_tree_file =          '%s_species_tree.newick'    % output_prefix

################################################### get species tree ###################################################

# get species tree
# Tests for presence of the tmp folder and deletes it
if os.path.exists(tmp_folder):
    os.system('rm -r ' + tmp_folder)
os.mkdir(tmp_folder)

# List all prokka dirs in the target folder
prokka_files = [i for i in os.listdir(path_to_prokka) if os.path.isdir(path_to_prokka + '/' + i)]
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Detected %i input genomes' % len(prokka_files))

# Running hmmsearch on each file
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running hmmsearch...')
for f in prokka_files:
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
                if (float(splitLine[13]) > hmm_score):
                    hmm_name = splitLine[0]
                    hmm_pos1 = int(splitLine[17]) - 1
                    hmm_pos2 = int(splitLine[18])
                    hmm_score = float(splitLine[13])
            else:
                file_out = open(tmp_folder + '/' + hmm_id + '.fasta', 'a+')
                file_out.write('>' + f + '\n')
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
            file_out = open(tmp_folder + '/' + hmm_id + '.fasta', 'a+')
            file_out.write('>' + f + '\n')
            if hmm_name != '':
                seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
            file_out.write(str(seq) + '\n')
            file_out.close()

# Call mafft to align all single fasta files with hmms
files = os.listdir(tmp_folder)
fastaFiles = [i for i in files if i.endswith('.fasta')]
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running mafft...')
for f in fastaFiles:
    fastaFile1 = '%s/%s' % (tmp_folder, f)
    fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
    os.system(
        pwd_mafft_exe + ' --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)

# concatenating the single alignments
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Concatenating alignments...')
concatAlignment = {}
for element in prokka_files:
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

    for element in prokka_files:
        if element in proteinSequence.keys():
            concatAlignment[element] += proteinSequence[element]
        else:
            concatAlignment[element] += '-' * alignmentLength

# writing alignment to file
file_out = open(combined_alignment_file, 'w')
for element in prokka_files:
    file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
file_out.close()

# calling fasttree for tree calculation
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running fasttree...')
os.system('%s -quiet %s > %s' % (pwd_fasttree_exe, combined_alignment_file, newick_tree_file))

# Decomment the two following lines if tree is rooted but should be unrooted
# phyloTree = dendropy.Tree.get(path='phylogenticTree.phy', schema='newick', rooting='force-unrooted')
# dendropy.Tree.write_to_path(phyloTree, 'phylogenticTree_unrooted.phy', 'newick')
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' The built species tree was exported to %s' % newick_tree_file)


######################################################## grouping ######################################################

# get bin to taxon dict
bin_to_taxon_dict = {}
if taxon_classification_file != None:
    bin_to_taxon_dict = get_rank_assignment_dict(selected_rank, taxon_classification_file)

# read in tree
tree_in = Tree(newick_tree_file, format=3)

# get sorted leaf node list
leaf_node_list = []
for leaf_node in tree_in:
    leaf_node_list.append(leaf_node.name)
leaf_node_list = sorted(leaf_node_list)

# get distance matrix from tree file
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get distance matrix from tree')

all_distances_lol = get_distance_matrix(newick_tree_file)

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
pwd_png_file_group =        '%s/%s'                     % (wd, png_file_group)
pwd_grouping_file =         '%s/%s'                     % (wd, grouping_file)
pwd_grouping_file_temp =    '%s/%s'                     % (wd, grouping_file_temp)

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
add_group_to_tree_R_cmd = 'Rscript %s -t %s -g %s > /dev/null' % (add_group_to_tree_R, newick_tree_file, pwd_grouping_file)
os.system(add_group_to_tree_R_cmd)

# report done
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping step done!')
