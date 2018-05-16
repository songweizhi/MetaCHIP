import os
import re
import glob
import shutil
import argparse
from sys import stdout
import configparser
from Bio import SeqIO
from ete3 import NCBITaxa, Tree
from PIL import Image, ImageDraw, ImageFont


class BinRecord(object):

    def __init__(self, name, id, group, group_without_underscore):
        self.name = name
        self.id = id
        self.group = group
        self.group_without_underscore = group_without_underscore


def plot_gene_tree(tree, tree_type, gene_name, tree_file_name, name_list, tree_image_folder):
    # set tree parameters
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    #ts.scale = 50
    ts.show_leaf_name = False
    tree_title = '%s (%s)' % (tree_type, gene_name)  # define tree title
    # tree title text setting
    ts.title.add_face(TextFace(tree_title,
                               fsize = 8,
                               fgcolor = 'black',
                               ftype = 'Arial',
                               tight_text = False),
                      column = 0)

    ts.rotation = 0  # from 0 to 360
    ts.show_scale = False
    ts.margin_top = 10  # top tree image margin
    ts.margin_bottom = 10  # bottom tree image margin
    ts.margin_left = 10  # left tree image margin
    ts.margin_right = 10  # right tree image margin
    ts.show_border = False  # set tree image border
    ts.branch_vertical_margin = 3  # 3 pixels between adjancent branches

    # set tree node style
    for each_node in tree.traverse():
        if each_node.is_leaf():  # leaf node parameters
            ns = NodeStyle()
            ns["shape"] = "circle"  # dot shape: circle, square or sphere
            ns["size"] = 0  # dot size
            ns['hz_line_width'] = 0.5  # branch line width
            ns['vt_line_width'] = 0.5  # branch line width
            ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
            ns['vt_line_type'] = 0  # branch line type
            if each_node.name in name_list:
                ns["fgcolor"] = "red"  # the dot setting
                # the node name text setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'red',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')
                each_node.set_style(ns)
            else:
                ns["fgcolor"] = "blue"  # the dot setting
                # the node name text setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'black',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')
                each_node.set_style(ns)
        else:  # non-leaf node parameters
            nlns = NodeStyle()
            nlns["size"] = 0  # dot size
            each_node.set_style(nlns)
    # set figures size
    tree.render('%s/%s_%s.png' % (tree_image_folder, tree_type, tree_file_name), w = 900, units = "px", tree_style = ts)


def plot_species_tree(tree_newick, tree_type, gene_name, tree_file_name, name_list, tree_image_folder):
    # set tree parameters
    tree = Tree(tree_newick, format = 8)
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    ts.show_leaf_name = False
    tree_title = tree_type + ' (' + gene_name + ')'  # define tree title
    # set tree title text parameters
    ts.title.add_face(TextFace(tree_title,
                               fsize = 8,
                               fgcolor = 'black',
                               ftype = 'Arial',
                               tight_text = False),
                      column = 0)  # tree title text setting
    # set layout parameters
    ts.rotation = 0  # from 0 to 360
    ts.show_scale = False
    ts.margin_top = 10  # top tree image margin
    ts.margin_bottom = 10  # bottom tree image margin
    ts.margin_left = 10  # left tree image margin
    ts.margin_right = 10  # right tree image margin
    ts.show_border = False  # set tree image border
    ts.branch_vertical_margin = 3  # 3 pixels between adjancent branches

    # set tree node style
    for each_node in tree.traverse():
        # leaf node parameters
        if each_node.is_leaf():
            ns = NodeStyle()
            ns['shape'] = 'circle'  # dot shape: circle, square or sphere
            ns['size'] = 0  # dot size
            ns['hz_line_width'] = 0.5  # branch line width
            ns['vt_line_width'] = 0.5  # branch line width
            ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
            ns['vt_line_type'] = 0  # branch line type
            if each_node.name in name_list:
                ns['fgcolor'] = 'red'  # the dot setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'red',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')  # the node name text setting
                each_node.set_style(ns)
            else:
                ns['fgcolor'] = 'blue'  # the dot setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize = 8,
                                            fgcolor = 'black',
                                            tight_text = False,
                                            bold = False),
                                   column = 0,
                                   position = 'branch-right')  # the node name text setting
                each_node.set_style(ns)

        # non-leaf node parameters
        else:
            nlns = NodeStyle()
            nlns['size'] = 0  # dot size
            each_node.add_face(TextFace(each_node.name,
                                        fsize = 4,
                                        fgcolor = 'black',
                                        tight_text = False,
                                        bold = False),
                               column = 5,
                               position = 'branch-top')  # non-leaf node name text setting)
            each_node.set_style(nlns)
    # set figures size
    tree.render('%s/%s_%s.png' % (tree_image_folder, tree_type, tree_file_name), w = 900, units = 'px', tree_style = ts)



############################################## Read in configuration file ##############################################

parser = argparse.ArgumentParser()
config = configparser.ConfigParser()

parser.add_argument('-cfg',
                    required=True,
                    help='path to configuration file')

args = vars(parser.parse_args())
pwd_cfg_file = args['cfg']
config.read(pwd_cfg_file)

ffn_file = config['FILES_AND_PARAMETERS']['ffn_file']
grouping_file = config['FILES_AND_PARAMETERS']['grouping_file']
cover_cutoff = config['FILES_AND_PARAMETERS']['cover_cutoff']
align_len_cutoff = config['FILES_AND_PARAMETERS']['align_len_cutoff']
identity_percentile = int(config['FILES_AND_PARAMETERS']['identity_percentile'])
ending_match_length = int(config['FILES_AND_PARAMETERS']['ending_match_length'])
ortholog_group_folder_name = config['FILES_AND_PARAMETERS']['ortholog_group_folder_name']

pwd_phyml_exe = config['DEPENDENCIES']['path_to_phyml_executable']
pwd_muscle_exe = config['DEPENDENCIES']['path_to_muscle_executable']
pwd_gblocks_exe = config['DEPENDENCIES']['path_to_gblocks_executable']
pwd_clustalo_exe = config['DEPENDENCIES']['path_to_clustalo_executable']
pwd_fasttree_exe = config['DEPENDENCIES']['path_to_fasttree_executable']
pwd_alignment_filter_script = config['DEPENDENCIES']['path_to_alignment_filter_script']
programs_for_HGT_prediction = config['DEPENDENCIES']['programs_for_HGT_prediction'].split(' ')
pwd_AnGST_exe = config['DEPENDENCIES']['path_to_AnGST_executable']
pwd_ranger_exe = config['DEPENDENCIES']['path_to_Ranger_executable']

############################################### Define folder/file name ################################################

wd = os.getcwd()
op_folder = 'output_ip' + str(identity_percentile) + '%_al' + str(align_len_cutoff) + 'bp_c' + str(cover_cutoff) + '%_e' + str(ending_match_length) + 'bp'

gene_tree_seq_folder =       'gene_tree_seq'
gene_tree_folder_ranger =    'gene_tree_ranger_txt'
gene_tree_folder_angst =     'gene_tree_angst_txt'
species_tree_folder_plot =   'species_tree_plot'
species_tree_folder_angst =  'species_tree_angst'
species_tree_folder_ranger = 'species_tree_ranger'
ranger_inputs_folder_name =  'Ranger_input'
ranger_outputs_folder_name = 'Ranger_output'
candidates_file_name =       'HGT_candidates.txt'
ranger_wd_name =             'Ranger-DTL_wd'
angst_inputs_folder_name =   'AnGST_wd'
output_tree_folder_name =    'Tree'
tree_image_folder_name =     'combined_tree_images'
#qsub_folder_for_AnGST = 'qsub_files_for_AnGST'
#angst_output_folder_name = 'gene_tree_angst_outputs'

pwd_grouping_file =              '%s/%s'         % (wd, grouping_file)
pwd_ffn_file =                   '%s/%s'         % (wd, ffn_file)
pwd_ortholog_group_folder =      '%s/%s'         % (wd, ortholog_group_folder_name)
pwd_candidates_file =            '%s/%s/%s'      % (wd, op_folder, candidates_file_name)
pwd_ranger_wd =                  '%s/%s/%s/%s/'  % (wd, op_folder, output_tree_folder_name, ranger_wd_name)
pwd_angst_inputs_folder =        '%s/%s/%s/%s'   % (wd, op_folder, output_tree_folder_name, angst_inputs_folder_name)
pwd_ranger_inputs_folder =       '%s/%s'         % (pwd_ranger_wd, ranger_inputs_folder_name)
pwd_ranger_outputs_folder =      '%s/%s'         % (pwd_ranger_wd, ranger_outputs_folder_name)
pwd_op_tree_folder =             '%s/%s/%s'      % (wd, op_folder, output_tree_folder_name)
pwd_tree_image_folder =          '%s/%s'         % (pwd_op_tree_folder, tree_image_folder_name)
pwd_gene_tree_seq_folder =       '%s/%s'         % (pwd_op_tree_folder, gene_tree_seq_folder)
pwd_gene_tree_folder_ranger =    '%s/%s'         % (pwd_op_tree_folder, gene_tree_folder_ranger)
pwd_species_tree_folder_plot =   '%s/%s'         % (pwd_op_tree_folder, species_tree_folder_plot)
pwd_species_tree_folder_angst =  '%s/%s'         % (pwd_op_tree_folder, species_tree_folder_angst)
pwd_species_tree_folder_ranger = '%s/%s'         % (pwd_op_tree_folder, species_tree_folder_ranger)

########################################################################################################################


# forward to working directory
os.chdir(wd)


# get list of match pair list
candidates = open(pwd_candidates_file)
candidates_list = []
for match_group in candidates:
    match_group_split = match_group.strip().split('\t')
    #match_group_split = sorted(match_group_split)
    candidates_list.append(match_group_split)


# get all ortholog groups
clusters_original = [os.path.basename(file_name) for file_name in glob.glob('%s/*.fna' % pwd_ortholog_group_folder)]
# remove " ' " from ortholog group name
clusters = []
for cluster_o in clusters_original:
    if "\'" in cluster_o:
        cmd_line_name = cluster_o.replace('\'', '\\\'')
        cluster_new = cluster_o.replace('\'', '')
        clusters.append(cluster_new)
        os.system('mv %s/%s %s/%s' % (pwd_ortholog_group_folder, cmd_line_name, pwd_ortholog_group_folder, cluster_new))
    else:
        clusters.append(cluster_o)


# get dict to hold members of each ortholog group
clusters_dict = {}
for cluster in clusters:
    cluster_key = cluster.split('.')[0]
    cluster_value = []
    members = open(pwd_ortholog_group_folder + '/' + cluster)
    for member in members:
        if member.startswith('>'):
            member_name = member.split('|')[0][4:-1]
            cluster_value.append(member_name)
    clusters_dict[cluster_key] = cluster_value


# get bin_record_list and genome name list
bin_taxon_ids = open(pwd_grouping_file)
bin_record_list = []
genome_name_list = []
name_to_group_number_dict = {}
name_to_group_number_without_underscore_dict = {}
bin_group_list = []
bin_group_without_underscore_list = []
for each_bin in bin_taxon_ids:
    each_bin_split = each_bin.strip().split(',')
    bin_group = each_bin_split[0]
    bin_group_without_underscore = bin_group.split('_')[0] + bin_group.split('_')[1]
    bin_name = each_bin_split[1]
    bin_id = each_bin_split[2]
    bin_parent_id = each_bin_split[3]
    name_to_group_number_dict[bin_name] = bin_group
    name_to_group_number_without_underscore_dict[bin_name] = bin_group_without_underscore
    bin_record = BinRecord(bin_name, bin_id, bin_group, bin_group_without_underscore)
    bin_record_list.append(bin_record)
    genome_name_list.append(bin_name)
    bin_group_list.append(bin_group)
    bin_group_without_underscore_list.append(bin_group_without_underscore)


##################################### Create folders to hold species/gene trees ######################################

# create output_tree_folder, gene_tree_seq and gene_tree_txt folder
if not os.path.isdir(pwd_op_tree_folder):
    os.mkdir(pwd_op_tree_folder)
    os.mkdir(pwd_gene_tree_seq_folder)
    os.mkdir(pwd_gene_tree_folder_ranger)
else:
    shutil.rmtree(pwd_op_tree_folder)
    os.mkdir(pwd_op_tree_folder)
    os.mkdir(pwd_gene_tree_seq_folder)
    os.mkdir(pwd_gene_tree_folder_ranger)

if 'Ranger-DTL' in programs_for_HGT_prediction:
    os.mkdir(pwd_species_tree_folder_ranger)
    os.mkdir(pwd_species_tree_folder_plot)

# if 'AnGST' in programs_for_HGT_prediction:
#     if not os.path.isdir(pwd_op_tree_folder + '/' + species_tree_folder_angst):
#         os.mkdir(pwd_op_tree_folder + '/' + species_tree_folder_angst)
#     else:
#         pass

####################################################### Main Code ######################################################

n = 1
for each_candidates in candidates_list:
    process_name = '___'.join(each_candidates)
    stdout.write("\rProcessing %sth of %s candidates: %s\n" % (str(n), str(len(candidates_list)), process_name))

    # get ortholog_list for each match pairs
    ortholog_list = []
    for each in clusters_dict:
        if (each_candidates[0] in clusters_dict[each]) or (each_candidates[1] in clusters_dict[each]):
            ortholog_list += clusters_dict[each]

    print(len(ortholog_list))
    if len(ortholog_list) == 0:
        print('No orthologes for the current HGT candidate')
    else:

        # get bin name list from candidates gene list
        each_candidates_bin_name = []
        for each_g in each_candidates:
            each_g_new = '_'.join(each_g.split('_')[:-1])
            each_candidates_bin_name.append(each_g_new)

        # uniq ortholog_list
        gene_member = []
        for each_g in ortholog_list:
            if each_g not in gene_member:
                gene_member.append(each_g)
        # get sequences of othorlog group to build gene tree
        seq_file_name = 'gene_tree_' + '___'.join(each_candidates) + '.seq'
        output_handle = open(pwd_gene_tree_seq_folder + '/' + seq_file_name, "w")
        records = SeqIO.parse(pwd_ffn_file, 'fasta')
        for record in records:
            if record.id in gene_member:
                record_id_split = record.id.split('_')
                bin_name_in_record_dot_name = '_'.join(record_id_split[:-1])
                for each_bin_record in bin_record_list:
                    if each_bin_record.name == bin_name_in_record_dot_name:
                        record.id = each_bin_record.group_without_underscore + '_' + record_id_split[-1]
                        name_len = len(record.id)
                        add_upp_num = 10 - name_len
                        record.id += add_upp_num * '_'
                        record.description = ''
                        SeqIO.write(record, output_handle, 'fasta')
        output_handle.close()


    ##################################################### tree builder #####################################################

        species_tree_ranger_file_name = 'species_tree_ranger_' + '___'.join(each_candidates)
        #species_tree_angst_file_name = 'species_tree_angst_' + '___'.join(each_candidates)
        gene_tree_file_name = 'gene_tree_' + '___'.join(each_candidates)

        # build tree with Ranger-DTL pipeline
        if 'Ranger-DTL' in programs_for_HGT_prediction:
            # prepare file name
            seq_file_name_non_extention = seq_file_name[:-4]
            first_clustalo_output = '%s.aln' % seq_file_name_non_extention
            filter_ranger_input_name = '%s.aln-gb' % seq_file_name_non_extention
            filter_ranger_output_name = '%s_aln_gb_filted.fasta' % seq_file_name_non_extention
            second_clustalo_output = '%s_aln_gb_filted.aln' % seq_file_name_non_extention
            fastTree_output = '%s.txt' % seq_file_name_non_extention

            # program parameters
            gblocks_parameters = '-t=d -b3=24 -b4=6 -b5=a'
            fasttree_parameters = '-gtr -nt'

            # path to input/output files
            path_to_first_clustalo_input = '%s/%s' % (pwd_gene_tree_seq_folder, seq_file_name)
            path_to_first_clustalo_output = '%s/%s' % (pwd_gene_tree_seq_folder, first_clustalo_output)
            path_to_ranger_filter_input = '%s/%s' % (pwd_gene_tree_seq_folder, filter_ranger_input_name)
            path_to_ranger_filter_output = '%s/%s' % (pwd_gene_tree_seq_folder, filter_ranger_output_name)
            path_to_second_clustalo_output = '%s/%s' % (pwd_gene_tree_seq_folder, second_clustalo_output)
            path_to_fasttree_output = '%s/%s' % (pwd_gene_tree_folder_ranger, fastTree_output)

            # prepare program command
            first_clustalo_cmd = '%s --force -i %s -o %s' % (pwd_clustalo_exe, path_to_first_clustalo_input, path_to_first_clustalo_output)
            gblocks_ranger_cmd = '%s %s %s' % (pwd_gblocks_exe, path_to_first_clustalo_output, gblocks_parameters)
            ranger_filter_cmd = 'python %s %s %s' % (pwd_alignment_filter_script, path_to_ranger_filter_input, path_to_ranger_filter_output)
            second_clustalo_cmd = '%s --force -i %s -o %s' % (pwd_clustalo_exe, path_to_ranger_filter_output, path_to_second_clustalo_output)
            FastTree_cmd = '%s %s %s > %s' % (pwd_fasttree_exe, fasttree_parameters, path_to_second_clustalo_output, path_to_fasttree_output)

            # execute commands
            print('Running first round clustal-o')
            os.system(first_clustalo_cmd)
            print('Running Gblocks')
            os.system(gblocks_ranger_cmd)
            print('Filtering alignment sequences')
            os.system(ranger_filter_cmd)
            print('Running second round clustal-o')
            os.system(second_clustalo_cmd)
            print('Running FastTree')
            os.system(FastTree_cmd)
            print('Got gene tree')

#     # # build tree with AnGST pipeline
#     # if 'AnGST' in programs_for_HGT_prediction:
#     #     # prepare file name
#     #     seq_file_name_non_extention = seq_file_name[:-4]
#     #     first_muscle_output_name = '%s_aln.fasta' % seq_file_name_non_extention
#     #     angst_filter_input_name = '%s_aln.fasta-gb' % seq_file_name_non_extention
#     #     angst_filter_output_name = '%s_aln_gb_filted.fasta' % seq_file_name_non_extention
#     #     second_muscle_output_name = '%s_aln_gb_filted_aln.phy' % seq_file_name_non_extention
#     #     phyml_output_name = '%s_aln_gb_filted_aln.phy_phyml_tree.txt' % seq_file_name_non_extention
#     #     phyml_tree_file_name = '%s.txt' % seq_file_name_non_extention
#     #
#     #     # program parameters
#     #     gblocks_parameters = '-t=d -b3=24 -b4=6 -b5=a'
#     #     phyml_parameters = '-a 1.0 -b 100 -o tl'
#     #
#     #     # path to input/output files
#     #     path_to_first_muscle_input_mac ='%s%s/%s' % (path_to_output_tree_folder, gene_tree_seq_folder, seq_file_name)
#     #     path_to_first_muscle_output_mac = '%s%s/%s' % (path_to_output_tree_folder, angst_output_folder_name, first_muscle_output_name)
#     #     path_to_angst_filter_input_mac = '%s%s/%s' % (path_to_output_tree_folder, angst_output_folder_name, angst_filter_input_name)
#     #     path_to_angst_filter_output_mac = '%s%s/%s' % (path_to_output_tree_folder, angst_output_folder_name, angst_filter_output_name)
#     #     path_to_second_muscle_output_mac = '%s%s/%s' % (path_to_output_tree_folder, angst_output_folder_name, second_muscle_output_name)
#     #     path_to_phyml_output_name_mac = '%s%s/%s' % (path_to_output_tree_folder, angst_output_folder_name, phyml_output_name)
#     #
#     #     # prepare program command
#     #     first_muscle_cmd = '%s -in %s -out %s' % (path_to_muscle_mac, path_to_first_muscle_input_mac, path_to_first_muscle_output_mac)
#     #     gblocks_angst_cmd = '%s %s %s' % (path_to_gblocks_mac, path_to_first_muscle_output_mac, gblocks_parameters)
#     #     angst_filter_cmd = 'python %s %s %s' % (path_to_alignment_filter_script_mac, path_to_angst_filter_input_mac, path_to_angst_filter_output_mac)
#     #     second_muscle_cmd = '%s -phyi -in %s -out %s' % (path_to_muscle_mac, path_to_angst_filter_output_mac, path_to_second_muscle_output_mac)
#     #     phyml_cmd = '%s %s -i %s' % (path_to_phyml_mac, phyml_parameters, path_to_second_muscle_output_mac)
#     #
#     #     # execute commands
#     #     os.system(first_muscle_cmd)
#     #     os.system(gblocks_angst_cmd)
#     #     os.system(angst_filter_cmd)
#     #     os.system(second_muscle_cmd)
#     #     os.system(phyml_cmd)

#################################################### Get species tree ##################################################

        # get bin name
        genome_member = []
        for each_gene in gene_member:
            # need to be changed
            each_gene_split = each_gene.split('_')
            genome_name = '_'.join(each_gene_split[:-1])
            genome_member.append(genome_name)
        print(gene_member)
        print(genome_member)

        # get  taxon id list for a sub-tree
        id_list = []
        for genome in genome_member:
            for each_bin_record in bin_record_list:
                if each_bin_record.name == genome:
                    id_list.append(each_bin_record.id)
        print(id_list)

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

        # change leaf name from bin_id to bin_name and write to output file
        tree = Tree(tree_phylo.write(format = 8))
        for leaf in tree.iter_leaves():
            for each_bin_record_2 in bin_record_list:
                if each_bin_record_2.id == leaf.name:
                    leaf.name = each_bin_record_2.group_without_underscore + '_' + each_bin_record_2.name

        if 'Ranger-DTL' in programs_for_HGT_prediction:
            species_tree_ranger_out = open(pwd_species_tree_folder_ranger + '/' + species_tree_ranger_file_name + '.txt', 'w')
            species_tree_ranger_out.write(tree.write(format = 8))
            species_tree_ranger_out.close()

        # if 'AnGST' in programs_for_HGT_prediction:
        #     species_tree_angst_out = open(pwd_op_tree_folder + '/' + species_tree_folder_angst + '/' + species_tree_angst_file_name + '.txt', 'w')
        #     # add root and its branch length to the tree, to fufill AnGST's requirement
        #     modified_species_tree = '(' + tree.write()[:-1] + ':1);'
        #     species_tree_angst_out.write(modified_species_tree)
        #     species_tree_angst_out.close()
#
    n += 1


########################################################################################################################
########################################################################################################################


# create folder to hold combined species tree and gene tree
if not os.path.exists(pwd_tree_image_folder):
    os.makedirs(pwd_tree_image_folder)
else:
    shutil.rmtree(pwd_tree_image_folder)
    os.makedirs(pwd_tree_image_folder)


# get members of all clusters
clusters = [os.path.basename(file_name) for file_name in glob.glob(pwd_ortholog_group_folder + '/*.fna')]
clusters_dict = {}
for cluster in clusters:
    cluster_key = cluster.split('.')[0]
    cluster_value = []
    members = open('%s/%s' % (pwd_ortholog_group_folder, cluster))
    for member in members:
        if member.startswith('>'):
            member_name = member.split('|')[0][4:-1]
            cluster_value.append(member_name)
    clusters_dict[cluster_key] = cluster_value


#################################### Prepare AnGST and Ranger-DTL working directory ####################################

# # prepare AnGST working directory
# if not os.path.exists(pwd_angst_inputs_folder):
#     os.makedirs(pwd_angst_inputs_folder)
# else:
#     shutil.rmtree(pwd_angst_inputs_folder)
#     os.makedirs(pwd_angst_inputs_folder)

# prepare Ranger-DTL working directory
if not os.path.exists(pwd_ranger_wd):
    os.makedirs(pwd_ranger_wd)
    os.makedirs(pwd_ranger_inputs_folder)
    os.makedirs(pwd_ranger_outputs_folder)
else:
    shutil.rmtree(pwd_ranger_wd)
    os.makedirs(pwd_ranger_wd)
    os.makedirs(pwd_ranger_inputs_folder)
    os.makedirs(pwd_ranger_outputs_folder)


####################################################### Main Code ######################################################

#AnGST_commands_file = open('%sAnGST_commands.txt' % pwd_angst_inputs_folder, 'w')

candidates_file_open = open(pwd_candidates_file)
candidates = []
for each_can in candidates_file_open:
    each_can_split = each_can.strip().split('\t')
    each_can_split = sorted(each_can_split)
    new = '%s___%s' % (each_can_split[0], each_can_split[1])
    candidates.append(new)
candidates = sorted(candidates)

processing_num = 1
for candidate in candidates:
    stdout.write("\rProcessing %sth of %s candidates: %s" % (str(processing_num), str(len(candidates)), candidate))

    # import species tree and gene tree
    ranger_species_tree = Tree('%s/species_tree_ranger_%s.txt' % (pwd_species_tree_folder_ranger, candidate), format = 1)
    ranger_gene_tree = Tree('%s/gene_tree_%s.txt' % (pwd_gene_tree_folder_ranger, candidate), format=1)
    species_tree_for_plot = '%s/species_tree_plot_%s.txt' % (pwd_species_tree_folder_plot, candidate)
    # assign new leaf name to gene tree
    for each_g_leaf in ranger_gene_tree:
        each_g_leaf_name_split = each_g_leaf.name.split('_')
        group_number = each_g_leaf_name_split[0]
        gene_id = each_g_leaf_name_split[1]
        for each_bin_record in bin_record_list:
            if group_number == each_bin_record.group_without_underscore:
                new_g_leaf_name = '%s_%s_%s' % (group_number, each_bin_record.name, gene_id)
                each_g_leaf.name = new_g_leaf_name

    # get leaf name style list
    candidate_split_gene = candidate.split('___')
    gene_name_with_group = []
    bin_name_with_group = []
    for each_candidate in candidate_split_gene:
        each_candidate_split = each_candidate.split('_')
        bin_name = '%s_%s' % (each_candidate_split[0], each_candidate_split[1])
        for bin_record in bin_record_list:
            if bin_record.name == bin_name:
                bin_name_group = '%s_%s' % (bin_record.group_without_underscore, bin_name)
                bin_name_with_group.append(bin_name_group)
                gene_name_group = '%s_%s' % (bin_record.group_without_underscore, each_candidate)
                gene_name_with_group.append(gene_name_group)

    # get gene name
    for each_cluster in clusters_dict:
        if candidate.split('___')[0] in clusters_dict[each_cluster]:
            gene_name = each_cluster
            # plot species tree and gene tree
            plot_species_tree(species_tree_for_plot, 'Species_Tree', gene_name, candidate, bin_name_with_group, pwd_tree_image_folder)
            plot_gene_tree(ranger_gene_tree, 'Gene_Tree', gene_name, candidate, gene_name_with_group, pwd_tree_image_folder)


#################################################### Run Ranger-DTL ####################################################

    # define Ranger-DTL input file name
    ranger_inputs_file_name = candidate + '.txt'
    ranger_outputs_file_name = candidate + '_output.txt'
    pwd_ranger_inputs = '%s/%s' % (pwd_ranger_inputs_folder, ranger_inputs_file_name)
    pwd_ranger_outputs = '%s/%s' % (pwd_ranger_outputs_folder, ranger_outputs_file_name)
    ranger_inputs_file = open(pwd_ranger_inputs, 'w')

    # change species tree leaf name for Ranger-DTL
    species_tree_ranger = ranger_species_tree
    for each_sr_leaf in species_tree_ranger:
        each_sr_leaf.name = each_sr_leaf.name.split('_')[0]

    # change gene tree leaf name for Ranger-DTL
    gene_tree_ranger = ranger_gene_tree
    for each_gr_leaf in gene_tree_ranger:
        each_gr_leaf_name_split = each_gr_leaf.name.split('_')
        each_gr_leaf.name = each_gr_leaf.name.split('_')[0]

    ranger_inputs_file.write('%s\n[&U]%s\n' % (species_tree_ranger.write(format = 9), gene_tree_ranger.write(format = 9)))
    ranger_inputs_file.close()

    # run Ranger-DTL
    ranger_parameters = '-q -D 2 -T 3 -L 1'
    ranger_cmd = '%s %s -i %s -o %s' % (pwd_ranger_exe, ranger_parameters, pwd_ranger_inputs, pwd_ranger_outputs)
    os.system(ranger_cmd)

    # parse prediction result
    ranger_result = open(pwd_ranger_outputs)
    predicted_transfers = []
    for each_line in ranger_result:
        if 'Transfer' in each_line:
            if not each_line.startswith('The minimum reconciliation cost'):
                mapping = each_line.strip().split(':')[1].split(',')[1]
                recipient = each_line.strip().split(':')[1].split(',')[2]
                donor_p = mapping.split('-->')[1][1:]
                donor_p_group = re.search('[a-zA-Z]*', donor_p).group()
                recipient_p = recipient.split('-->')[1][1:]
                recipient_p_group = re.search('[a-zA-Z]*', recipient_p).group()
                predicted_transfer = donor_p + '-->' + recipient_p
                if (mapping.split('-->')[1][1:] in bin_group_without_underscore_list) and (recipient.split('-->')[1][1:] in bin_group_without_underscore_list) and (donor_p_group != recipient_p_group):
                    predicted_transfers.append(predicted_transfer)
            else:
                pass

    # get two possible transfer situation
    candidate_split_group = []
    for each_gene in candidate_split_gene:
        each_gene_bin_name = '_'.join(each_gene.split('_')[:-1])
        for each_bin_record in bin_record_list:
            if each_gene_bin_name == each_bin_record.name:
                each_gene_g = each_bin_record.group
                candidate_split_group.append(each_gene_g)
    possible_hgt_1 = candidate_split_group[0].split('_')[0] + candidate_split_group[0].split('_')[1] + '-->' + candidate_split_group[1].split('_')[0] + candidate_split_group[1].split('_')[1]
    possible_hgt_2 = candidate_split_group[1].split('_')[0] + candidate_split_group[1].split('_')[1] + '-->' + candidate_split_group[0].split('_')[0] + candidate_split_group[0].split('_')[1]
    possible_hgts = [possible_hgt_1, possible_hgt_2]


####################### Combine Species/Gene Tree Together and add Ranger-DTL prediction results #######################

    # read in species/gene tree image
    species_tree_image = Image.open('%s/Species_Tree_%s.png' % (pwd_tree_image_folder, candidate))
    gene_tree_image = Image.open('%s/Gene_Tree_%s.png' % (pwd_tree_image_folder, candidate))
    images = [species_tree_image, gene_tree_image]

    # get new width and height
    widths = []
    heights = []
    text_height = (len(predicted_transfers)//9 + 1) * 40 + 80
    for image_size in images:
        width = image_size.size[0]
        height = image_size.size[1]
        widths.append(width)
        heights.append(height)
    new_width = sum(widths) + 90  # 30 left, 30 right and 30 in the middle
    new_height = max(heights) + 30 + text_height  # 30 in the top, text_height for the text

    # create a new image
    new_image = Image.new('RGB', color = (255, 255, 255), size = (new_width, new_height))  # setup background and size
    x_starting = 30
    for image_paste in images:
        new_image.paste(image_paste, (x_starting, 30))
        x_starting += image_paste.size[0] + 50

    # add Ranger-DTL prediction results to the combined image
    add_text = ImageDraw.Draw(new_image)
    # define font
    font_arial = ImageFont.truetype("/Library/Fonts/Arial.ttf", 24)
    font_arial_bold = ImageFont.truetype("/Library/Fonts/Arial Bold.ttf", 24)
    font_arial_black = ImageFont.truetype("/Library/Fonts/Arial Black.ttf", 24)
    # predicted_transfers
    x_start = 30
    y_start_t = 30 + max(heights)
    y_start = 30 + max(heights) + 40
    add_text.text((x_start, y_start_t), 'Ranger-DTL-U predicted HGTs (%s): ' % len(predicted_transfers), (0, 0, 0), font = font_arial_bold)


    for hgt in predicted_transfers:
        hgt_split = hgt.split('-->')
        hgt_d = hgt_split[0]
        hgt_r = hgt_split[1]

        if x_start < 1800:
            if hgt in possible_hgts:
                add_text.text((x_start, y_start), hgt, (225, 0, 0), font = font_arial)
                x_start += 200
            else:
                if (hgt_d in bin_group_without_underscore_list) and (hgt_r in bin_group_without_underscore_list) and (hgt_d[0] != hgt_r[0]):
                    add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                    x_start += 200
                else:
                    add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                    x_start += 200
        elif x_start >= 1800:
            x_start = 30
            y_start += 40
            if hgt in possible_hgts:
                add_text.text((x_start, y_start), hgt, (225, 0, 0), font = font_arial)
                x_start += 200
            else:
                if (hgt_d in bin_group_without_underscore_list) and (hgt_r in bin_group_without_underscore_list) and (hgt_d[0] != hgt_r[0]):
                    add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                    x_start += 200
                else:
                    add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                    x_start += 200

    new_image.save('%s/Combined_Tree_%s.png' % (pwd_tree_image_folder, candidate))
    os.remove('%s/Species_Tree_%s.png' % (pwd_tree_image_folder, candidate))  # remove species tree
    os.remove('%s/Gene_Tree_%s.png' % (pwd_tree_image_folder, candidate))  # remove gene tree

# ############################################### Prepare Input for AnGST ################################################
#
#     # generate folders for outputs
#     if not os.path.exists('%s/%s' % (pwd_angst_inputs_folder, candidate)):
#         os.makedirs('%s/%s' % (pwd_angst_inputs_folder, candidate))
#     else:
#         shutil.rmtree('%s/%s' % (pwd_angst_inputs_folder, candidate))
#         os.makedirs('%s/%s' % (pwd_angst_inputs_folder, candidate))
#
#     # read in AnGST species/gene tree
#     angst_species_tree = Tree('%s/species_tree_angst_%s.txt' % (pwd_species_tree_folder_angst, candidate), format = 1)
#     angst_gene_tree = Tree('%s/gene_tree_%s.txt' % (pwd_gene_tree_folder_ranger, candidate), format = 1)
#
#     # change leaf name in AnGST species tree
#     for each_asl in angst_species_tree:
#         each_asl.name = each_asl.name.split('_')[0]
#
#     # change leaf name in AnGST gene tree
#     for each_agl in angst_gene_tree:
#         each_agl_name_split = each_agl.name.split('_')
#         each_agl.name = '%s_%s' % (each_agl_name_split[0], each_agl_name_split[1])
#
#     # get species/gene tree for AnGST input
#     angst_species_tree_new = open('%s/%s/species_tree_%s.txt' % (pwd_angst_inputs_folder, candidate, candidate), 'w')
#     angst_gene_tree_new = open('%s/%s/gene_tree_%s.txt' % (pwd_angst_inputs_folder, candidate, candidate), 'w')
#     angst_species_tree_new.write(angst_species_tree.write() + '\n')
#     angst_gene_tree_new.write(angst_gene_tree.write() + '\n')
#     angst_species_tree_new.close()
#     angst_gene_tree_new.close()
#
#     # get AnGST.input file
#     AnGST_dot_input_file = open('%s/%s/AnGST.input' % (pwd_angst_inputs_folder, candidate), 'w')
#     line_1 = 'species=%s/%s/species_tree_%s.txt' % (pwd_angst_inputs_folder, candidate, candidate)
#     line_2 = 'gene=%s/%s/gene_tree_%s.txt' % (pwd_angst_inputs_folder, candidate, candidate)
#     line_3 = 'output=%s/%s/%s_out' % (pwd_angst_inputs_folder, candidate, candidate)
#     line_4 = 'penalties=%s/%s/penalty.file' % (pwd_angst_inputs_folder, candidate)
#     AnGST_dot_input_file.write('%s\n%s\n%s\n%s\n' % (line_1, line_2, line_3, line_4))
#     AnGST_dot_input_file.close()
#
#     # get penalty.file file
#     penalty_dot_file = open('%s/%s/penalty.file' % (pwd_angst_inputs_folder, candidate), 'w')
#     penalty_dot_file_content = 'hgt : 3.0\ndup : 2.0\nlos : 1.0\nspc : 0.0\n'
#     penalty_dot_file.write(penalty_dot_file_content)
#     penalty_dot_file.close()
#
#     # prepare commands to run AnGST
#     AnGST_dot_input_file_location = '%s/%s/AnGST.input' % (pwd_angst_inputs_folder, candidate)
#     #AnGST_commands_file.write('python %s %s\n' % (path_to_AnGST_executable, AnGST_dot_input_file_location))
#     processing_num += 1
#
# # Close files
# #AnGST_commands_file.close()

print('\nAll done!')
