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


def get_node_distance(tree, node_1, node_2):
    distance = tree.get_distance(node_1, node_2)
    return distance


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


pwd_newick_tree_file = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05_species_tree.newick'
pwd_newick_tree_file = '/Users/songweizhi/Desktop/human_gut_species_tree.newick'

# read in tree
tree_in = Tree(pwd_newick_tree_file, format=0)
#tree_in = Tree(pwd_newick_tree_file, format=0)


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

max_d = None
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

print(group_number)









output_prefix = 'test_Tree'





# define output file name with grouping included
png_file_group =            '%s_grouping_g%s.png'       % (output_prefix, group_number)
grouping_file =             '%s_grouping_g%s.txt'       % (output_prefix, group_number)
grouping_file_temp =        '%s_grouping_g%s_tmp.file'  % (output_prefix, group_number)
pwd_png_file_group =        '%s/%s'                     % ('/Users/songweizhi/Desktop', png_file_group)
pwd_grouping_file =         '%s/%s'                     % ('/Users/songweizhi/Desktop', grouping_file)
pwd_grouping_file_temp =    '%s/%s'                     % ('/Users/songweizhi/Desktop', grouping_file_temp)

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

# # calculate full dendrogram
# sleep(0.5)
# print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get plot for visualization')
# if taxon_classification_file != None:
#     plot_clustering_dendrogram(cluster, leaf_font_size, get_taxon, max_d, pwd_png_file_group)
# else:
#     plot_clustering_dendrogram(cluster, leaf_font_size, get_group, max_d, pwd_png_file_group)
#
# # call R
# current_wd = os.getcwd()
# os.chdir('/Users/songweizhi/Desktop')
# add_group_to_tree_R_cmd = 'Rscript %s -t %s -g %s -l %s > /dev/null' % (add_group_to_tree_R, newick_tree_file, grouping_file, label_shift)
# print(add_group_to_tree_R_cmd)
# os.system(add_group_to_tree_R_cmd)
# os.chdir(current_wd)
#
# # report done
# sleep(0.5)
# print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping step done!')
#
# sleep(0.5)
# print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' You may want to modify the grouping results based on the taxonomy\nof your input genomes/bins, you can do this by changing their group assignment specified\nin the first column of %s' % grouping_file)
#
