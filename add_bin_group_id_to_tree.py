
import os
import argparse
import itertools
from ete3 import Tree
from string import ascii_uppercase


usage = """

python /srv/scratch/z5039045/Scripts/add_bin_group_id_to_tree.py -t species_tree.newick -g grouping.txt -o species_tree_with_grouping.newick

"""





def get_id2name_dict(grouping_file):
    id2name_dict = {}
    for each in open(grouping_file):
        each_split = each.strip().split(',')
        group_id = each_split[0]
        genome_name = each_split[1]
        id2name_dict[genome_name] = group_id
    return id2name_dict


def add_grouping_to_tree(tree_file, id2name_dict, output_file):
    # rename leaf node name
    tree_file_name, tree_file_ext = os.path.splitext(tree_file)
    gene_tree = Tree(tree_file, format=3)
    gene_tree.resolve_polytomy(recursive=True)  # solving multifurcations

    # change gene tree leaf name for Ranger-DTL
    for leaf_node in gene_tree:
        leaf_node_name_new = '%s_%s' % (id2name_dict[leaf_node.name], leaf_node.name)
        leaf_node.name = leaf_node_name_new

    # write tree to new file
    tree_file_with_grouping_handle = open(output_file, 'w')
    tree_file_with_grouping_handle.write('%s\n' % (gene_tree.write(format=3)))
    tree_file_with_grouping_handle.close()


parser = argparse.ArgumentParser()

parser.add_argument('-t',
                    required=True,
                    help='tree file')

parser.add_argument('-g',
                    required=True,
                    help='grouping file')

parser.add_argument('-o',
                    required=True,
                    help='output tree with grouping')

args = vars(parser.parse_args())

tree_file = args['t']
output_file = args['o']
grouping_file = args['g']


id2name_dict = get_id2name_dict(grouping_file)
add_grouping_to_tree(tree_file, id2name_dict, output_file)
