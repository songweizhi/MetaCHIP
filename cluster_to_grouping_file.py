#!/usr/bin/env python

import os
import argparse
import itertools
from ete3 import Tree
from string import ascii_uppercase


usage = """

python /srv/scratch/z5039045/Scripts/cluster_to_grouping_file.py -t species_tree.newick -c grouping_out.txt -g grouping.txt

"""

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


def cluster_2_grouping_file(cluster_file, grouping_file):

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
    for each in open(t1_sorted):
        cluster_name = each.strip().split(',')[0]
        genome_name = each.strip().split(',')[1]

        if current_cluster_name == '':
            current_cluster_name = cluster_name
            grouping_file_handle.write('%s,%s\n' % (group_index_list[group_index_no], genome_name))
        elif current_cluster_name == cluster_name:
            grouping_file_handle.write('%s,%s\n' % (group_index_list[group_index_no], genome_name))
        elif current_cluster_name != cluster_name:
            current_cluster_name = cluster_name
            group_index_no += 1
            grouping_file_handle.write('%s,%s\n' % (group_index_list[group_index_no], genome_name))

    os.remove(t1)
    os.remove(t1_sorted)


def get_id2name_dict(grouping_file):
    id2name_dict = {}
    for each in open(grouping_file):
        each_split = each.strip().split(',')
        group_id = each_split[0]
        genome_name = each_split[1]
        id2name_dict[genome_name] = group_id
    return id2name_dict


def add_grouping_to_tree(tree_file, id2name_dict):
    # rename leaf node name
    tree_file_name, tree_file_ext = os.path.splitext(tree_file)
    tree_file_with_grouping = '%s_with_grouping%s' % (tree_file_name, tree_file_ext)
    gene_tree = Tree(tree_file, format=3)
    gene_tree.resolve_polytomy(recursive=True)  # solving multifurcations

    # change gene tree leaf name for Ranger-DTL
    # for leaf_node in gene_tree:
    #     leaf_node_name_new = '%s_%s' % (id2name_dict[leaf_node.name], leaf_node.name)
    #     leaf_node.name = leaf_node_name_new

    for leaf_node in gene_tree:
        leaf_node_name_new = '%s' % (id2name_dict[leaf_node.name])
        leaf_node.name = leaf_node_name_new


    # write tree to new file
    tree_file_with_grouping_handle = open(tree_file_with_grouping, 'w')
    tree_file_with_grouping_handle.write('%s\n' % (gene_tree.write(format=3)))
    tree_file_with_grouping_handle.close()


parser = argparse.ArgumentParser()

parser.add_argument('-c',
                    required=True,
                    help='cluster file')

parser.add_argument('-t',
                    required=True,
                    help='tree file')

parser.add_argument('-g',
                    required=True,
                    help='grouping file')

args = vars(parser.parse_args())

tree_file = args['t']
cluster_file = args['c']
grouping_file = args['g']


cluster_2_grouping_file(cluster_file, grouping_file)
id2name_dict = get_id2name_dict(grouping_file)
add_grouping_to_tree(tree_file, id2name_dict)
