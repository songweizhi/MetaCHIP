import os
import copy
import argparse
from Bio import Phylo
from datetime import datetime


def check_to_keep(clade_name_str, identified_taxon_list):
    # remove colon form clade name
    clade_name_str_no_colon = clade_name_str
    if ':' in clade_name_str_no_colon:
        clade_name_str_no_colon = clade_name_str_no_colon.split(':')[1]

    # split clade name if there are ';'
    for_check = []
    if ';' in clade_name_str_no_colon:
        clade_name_str_no_colon_split = clade_name_str_no_colon.split(';')

        # remove space at the begining or end for each split
        clade_name_str_no_colon_split_no_space = []
        for clade_name in clade_name_str_no_colon_split:
            if clade_name[0] == ' ':
                clade_name = clade_name[1:]
            if clade_name[-1] == ' ':
                clade_name = clade_name[:-1]
            clade_name_str_no_colon_split_no_space.append(clade_name)
        for_check = clade_name_str_no_colon_split_no_space
    else:
        for_check = [clade_name_str_no_colon]

    # check to keep or not
    clade_to_keep = 0
    for identified_taxon in identified_taxon_list:
        if identified_taxon in for_check:
            clade_to_keep = 1

    return clade_to_keep


def remove_unwanted_leaf_nodes(tree, identified_taxon_list):
    # copy tree
    tree_copy = copy.deepcopy(tree)

    removed_leaf_num = 0
    all_leaf_nodes = tree_copy.get_terminals()
    for leaf_node in all_leaf_nodes:
        leaf_node_name_str = str(leaf_node.name)
        leaf_node_to_keep = check_to_keep(leaf_node_name_str, identified_taxon_list)

        if leaf_node_to_keep == 0:
            tree_copy.collapse(leaf_node)
            removed_leaf_num += 1

    return tree_copy, removed_leaf_num


################################################# input #################################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-tree', dest='TREE', nargs='?', required=True,  type=str, help='GTDB tree file')
required.add_argument('-taxon', dest='TAXON', nargs='?', required=True, type=str, help='PrepIn produced group_to_taxon file')
required.add_argument('-output', dest='OUTPUT', nargs='?', required=True,  type=str, help='Output tree subset')

args = vars(parser.parse_args())
tree_file_in = args['TREE']
group_to_taxon_file = args['TAXON']
tree_file_out = args['OUTPUT']

# define tmp file name
tree_file_tmp_1 = '%s.tmp_1.tree' % tree_file_out
tree_file_tmp_2 = '%s.tmp_2.tree' % tree_file_out

time_format = '[%Y-%m-%d %H:%M:%S] '


################################################ store input information ###############################################

# read in tree
tree_in = Phylo.read(tree_file_in, 'newick')

# read in all identified taxons
identified_taxon_list = set()
for each_group in open(group_to_taxon_file):
    identified_taxon_list.add(each_group.strip().split(',')[1])
print(datetime.now().strftime(time_format) + 'The number of provided taxon: %s' % len(identified_taxon_list))


########################################## remove unwanted nodes recursively ###########################################

# remove unwanted nodes recursively
print(datetime.now().strftime(time_format) + 'Recursively removing unwanted nodes')
deleted_leaf_num = 1
n = 0
tree_in_copy = copy.deepcopy(tree_in)
while deleted_leaf_num > 0:
    tree_in_copy, deleted_leaf_num = remove_unwanted_leaf_nodes(tree_in_copy, identified_taxon_list)
    n += 1
    print(datetime.now().strftime(time_format) + 'Removed %s nodes in %sth round' % (deleted_leaf_num, n))

# write out tree
Phylo.write(tree_in_copy, tree_file_tmp_1, 'newick')


############################################# remove "100:" in clade name ##############################################

# read in tree
tree_tmp_1 = Phylo.read(tree_file_tmp_1, 'newick')
tree_tmp_1_copy = copy.deepcopy(tree_tmp_1)

for clade in tree_tmp_1_copy.find_clades():
    clade_name = str(clade.name)
    if ':' in clade_name:
        clade.name = clade_name.split(':')[1]

Phylo.write(tree_tmp_1_copy, tree_file_tmp_2, 'newick')


################################################ rename leaf nodes name ################################################

# read in tree
tree_tmp_2 = Phylo.read(tree_file_tmp_2, 'newick')
tree_tmp_2_copy = copy.deepcopy(tree_tmp_2)

# get all leaf nodes
all_leaf_nodes = tree_tmp_2_copy.get_terminals()
for leaf_node in all_leaf_nodes:
    leaf_node_name_str = str(leaf_node.name)

    if ';' in leaf_node_name_str:
        leaf_node_name_split = leaf_node_name_str.split(';')

        # remove space at the begining or end
        leaf_node_name_split_no_space = []
        for each_name in leaf_node_name_split:
            if each_name[0] == ' ':
                each_name = each_name[1:]
            if each_name[-1] == ' ':
                each_name = each_name[:-1]
            leaf_node_name_split_no_space.append(each_name)

        leaf_node_name_new = ''
        for identified_taxon in identified_taxon_list:
            if identified_taxon in leaf_node_name_split_no_space:
                leaf_node_name_new = identified_taxon

        leaf_node.name = leaf_node_name_new

# write out tree
Phylo.write(tree_tmp_2_copy, tree_file_out, 'newick')


# report
print(datetime.now().strftime(time_format) + 'Tree subset exported to: %s' % tree_file_out)


# print warning message if some provided node(s) were not found
extracted_leaf_nodes = tree_tmp_2_copy.get_terminals()
if len(extracted_leaf_nodes) < len(identified_taxon_list):

    extracted_leaf_node_list = []
    for extracted_leaf_node in extracted_leaf_nodes:
        extracted_leaf_node_list.append(str(extracted_leaf_node.name))

    un_extracted_nodes = []
    for provided_node in identified_taxon_list:
        if provided_node not in extracted_leaf_node_list:
            un_extracted_nodes.append(provided_node)

    print(datetime.now().strftime(time_format) + 'Warning!!! Found %s of %s provided nodes, missed: %s' % (len(extracted_leaf_nodes), len(identified_taxon_list), ', '.join(un_extracted_nodes)))


################################################### remove tmp files ###################################################

# remove tmp files
os.remove(tree_file_tmp_1)
os.remove(tree_file_tmp_2)
