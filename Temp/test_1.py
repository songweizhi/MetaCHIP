import os
import re
import glob
import shutil
from sys import stdout
import configparser
from Bio import SeqIO
from ete3 import NCBITaxa, Tree
from PIL import Image, ImageDraw, ImageFont
from lib.tree_ploter import plot_species_tree, plot_gene_tree
from lib.Bin_record_class import BinRecord


# get all id list for a sub-tree
id_list = []


# get species tree
ncbi = NCBITaxa()
tree_phylo = ncbi.get_topology(id_list, intermediate_nodes = 0)

# customize species tree for plot
tree_phylo_plot = ncbi.get_topology(id_list, intermediate_nodes = 1)
for each_node in tree_phylo_plot.traverse():
    node_name_list = []
    node_name_list.append(each_node.name)
    if node_name_list == ['']:
        pass
    else:
        if each_node.is_leaf():
            # change bin id to bin name
            each_node.name = ncbi.get_taxid_translator(node_name_list)[int(each_node.name)]
            # add group information to bin name
            each_node.name = name_to_group_number_without_underscore_dict[each_node.name] + '_' + each_node.name

        # for non-leaf node name, only display the first 12 characters
        else:
            if len(ncbi.get_taxid_translator(node_name_list)[int(each_node.name)]) <= 12:
                each_node.name = ncbi.get_taxid_translator(node_name_list)[int(each_node.name)]
            else:
                each_node.name = ncbi.get_taxid_translator(node_name_list)[int(each_node.name)][0:12] + '.'

species_tree_plot_out = open(pwd_species_tree_folder_plot + '/species_tree_plot_' + process_name + '.txt', 'w')
species_tree_plot_out.write(tree_phylo_plot.write(format = 8) + '\n')
species_tree_plot_out.close()