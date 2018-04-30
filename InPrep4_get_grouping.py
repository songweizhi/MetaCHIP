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
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import dendrogram


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

parser.add_argument('-t',
                    required=True,
                    help='tree file')

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

args = vars(parser.parse_args())
tree_file = args['t']
max_d = args['dc']
leaf_font_size = args['fs']
taxon_classification_file = args['taxon']
selected_rank = 'c'

########################################################################################################################

# disable warnings
warnings.filterwarnings("ignore")

# get bin to taxon dict
bin_to_taxon_dict = {}
if taxon_classification_file != None:
    bin_to_taxon_dict = get_rank_assignment_dict(selected_rank, taxon_classification_file)

# read in tree
tree_in = Tree(tree_file, format=3)

# get sorted leaf node list
leaf_node_list = []
for leaf_node in tree_in:
    leaf_node_list.append(leaf_node.name)
leaf_node_list = sorted(leaf_node_list)

# get distance matrix from tree file
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get distance matrix from tree')

all_distances_lol = get_distance_matrix(tree_file)

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

# define output file name
# input files
wd = os.getcwd()

tree_file_path, tree_file_name = os.path.split(tree_file)
tree_file_name_no_extension, tree_file_extension = os.path.splitext(tree_file_name)

png_file_group = '%s_grouping_g%s.png' % (tree_file_name_no_extension, group_number)
grouping_file = '%s_grouping_g%s.txt' % (tree_file_name_no_extension, group_number)
grouping_file_temp = '%s_grouping_g%s_tmp.file' % (tree_file_name_no_extension, group_number)

pwd_png_file_group = '%s/%s' % (wd, png_file_group)

pwd_grouping_file = '%s/%s' % (wd, grouping_file)
pwd_grouping_file_temp = '%s/%s' % (wd, grouping_file_temp)


# get grouping file
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Write out grouping file')
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

# calculate full dendrogram
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get plot for visualization')
if taxon_classification_file != None:
    plot_clustering_dendrogram(cluster, leaf_font_size, get_taxon, max_d, pwd_png_file_group)
else:
    plot_clustering_dendrogram(cluster, leaf_font_size, get_group, max_d, pwd_png_file_group)

# report done
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping finished')

add_group_to_tree_R = '~/R_scripts/newick_tree/add_group_to_tree.R'

add_group_to_tree_R_cmd = 'Rscript %s -t %s -g %s > /dev/null' % (add_group_to_tree_R, tree_file, pwd_grouping_file)

os.system(add_group_to_tree_R_cmd)

