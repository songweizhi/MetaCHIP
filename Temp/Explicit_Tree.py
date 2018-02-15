#!/usr/bin/env python
import os
import re
import glob
import shutil
import argparse
from ConfigParser import SafeConfigParser
from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


# to-do
# what if hiden files are more than one (section: get species tree for all input genomes)
# ortholog_list == gene_member, need to modify
# add gene name instead of concatenated candidates id to the species tree and gene tree
# what if not plot, how to prepare outputs
# what does checkm tree do
# for orthology groups, file name should be the protein name
# no ''' in input file name
# plot inter node of species tree


class BinRecord(object):

    def __init__(self, name, group, group_without_underscore):
        self.name = name
        self.group = group
        self.group_without_underscore = group_without_underscore


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def get_species_tree_alignment(tmp_folder, path_to_prokka, path_to_hmm, pwd_hmmsearch_exe, pwd_mafft_exe):

    # Tests for presence of the tmp folder and deletes it
    if os.path.exists(tmp_folder):
        os.system('rm -r ' + tmp_folder)
    os.mkdir(tmp_folder)

    # List all prokka dirs in the target folder
    prokka_files = [i for i in os.listdir(path_to_prokka) if os.path.isdir(path_to_prokka + '/' + i)]
    print('Detected %i input genomes' % len(prokka_files))

    # Running hmmsearch on each file
    print('Running hmmsearch...')
    for f in prokka_files:
        # call hmmsearch
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
                    if(float(splitLine[13]) > hmm_score):
                        hmm_name = splitLine[0]
                        hmm_pos1 = int(splitLine[17]) - 1
                        hmm_pos2 = int(splitLine[18])
                        hmm_score = float(splitLine[13])
                else:
                    file_out = open(tmp_folder + '/' + hmm_id + '.fasta', 'a+')
                    file_out.write('>' + f + '\n')
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
                seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
                file_out.write(str(seq) + '\n')
                file_out.close()

    # Call mafft to align all single fasta files with hmms
    files = os.listdir(tmp_folder)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    print('Running mafft...')
    for f in fastaFiles:
        fastaFile1 = '%s/%s' % (tmp_folder, f)
        fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
        os.system(pwd_mafft_exe + ' --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)


def get_species_tree_newick(tmp_folder, each_subset, pwd_fasttree_exe, tree_folder, tree_name):
    # concatenating the single alignments
    concatAlignment = {}
    for element in each_subset:
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

        for element in each_subset:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    pwd_alignment_file = '%s/%s.aln' % (tree_folder, tree_name)
    pwd_newick_file = '%s/%s.newick' % (tree_folder, tree_name)

    file_out = open(pwd_alignment_file, 'w')
    # for element in prokka_files:
    for element in each_subset:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()

    # calling fasttree for tree calculation
    os.system('%s -quiet %s > %s' % (pwd_fasttree_exe, pwd_alignment_file, pwd_newick_file))

    # Decomment the two following lines if tree is rooted but should be unrooted
    # phyloTree = dendropy.Tree.get(path='phylogenticTree.phy', schema='newick', rooting='force-unrooted')
    # dendropy.Tree.write_to_path(phyloTree, 'phylogenticTree_unrooted.phy', 'newick')
    # os.system('rm -r ' + tmp_folder)

    # draw species tree
    #species_tree = Phylo.read('phylogenticTree.phy', 'newick')
    #Phylo.draw_ascii(species_tree)


def plot_gene_tree(tree_newick, tree_type, gene_name, tree_file_name, name_list, tree_image_folder):

    # set tree parameters
    tree = Tree(tree_newick, format=1)
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    #ts.scale = 50
    ts.show_leaf_name = False
    tree_title = '%s (%s)' % (tree_type, gene_name)  # define tree title
    # tree title text setting
    ts.title.add_face(TextFace(tree_title,
                               fsize=8,
                               fgcolor='black',
                               ftype='Arial',
                               tight_text=False),
                      column=0)

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
                                            fsize=8,
                                            fgcolor='red',
                                            tight_text=False,
                                            bold=False),
                                   column=0,
                                   position='branch-right')
                each_node.set_style(ns)
            else:
                ns["fgcolor"] = "blue"  # the dot setting
                # the node name text setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize=8,
                                            fgcolor='black',
                                            tight_text=False,
                                            bold=False),
                                   column=0,
                                   position='branch-right')
                each_node.set_style(ns)
        else:  # non-leaf node parameters
            nlns = NodeStyle()
            nlns["size"] = 0  # dot size
            each_node.set_style(nlns)
    # set figures size
    tree.render('%s/%s.png' % (tree_image_folder, tree_file_name), w=900, units="px", tree_style=ts)


def plot_species_tree(tree_newick, tree_type, gene_name, tree_file_name, name_list, tree_image_folder):
    # set tree parameters
    tree = Tree(tree_newick, format=2)
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
    tree.render('%s/%s.png' % (tree_image_folder, tree_file_name), w = 900, units = 'px', tree_style = ts)


############################################## Read in configuration file ##############################################

parser = argparse.ArgumentParser()

parser.add_argument('-cfg',
                    required=True,
                    help='path to configuration file')

args = vars(parser.parse_args())
pwd_cfg_file = args['cfg']

# check whether config file exist
if not os.path.isfile(pwd_cfg_file):
    print('config file not found')
    exit()

config = SafeConfigParser()
config.read(pwd_cfg_file)

grouping_file =              config.get('FILES_AND_PARAMETERS', 'grouping_file')
cover_cutoff =               config.get('FILES_AND_PARAMETERS', 'cover_cutoff')
align_len_cutoff =           config.get('FILES_AND_PARAMETERS', 'align_len_cutoff')
identity_percentile =    int(config.get('FILES_AND_PARAMETERS', 'identity_percentile'))
ending_match_length =    int(config.get('FILES_AND_PARAMETERS', 'ending_match_length'))
prokka_output =              config.get('FILES_AND_PARAMETERS', 'prokka_outputs')
ortholog_group_folder_name = config.get('FILES_AND_PARAMETERS', 'orthologs_folder')
plot_tree =              int(config.get('FILES_AND_PARAMETERS', 'plot_tree'))
pwd_phylo_hmm =              config.get('FILES_AND_PARAMETERS', 'phylo_hmm')
pwd_ranger_exe =             config.get('DEPENDENCIES', 'path_to_ranger_executable')
pwd_hmmsearch_exe =          config.get('DEPENDENCIES', 'path_to_hmmsearch_executable')
pwd_mafft_exe =              config.get('DEPENDENCIES', 'path_to_mafft_executable')
pwd_fasttree_exe =           config.get('DEPENDENCIES', 'path_to_fasttree_executable')
# pwd_gblocks_exe = config['DEPENDENCIES']['path_to_gblocks_executable']
# pwd_alignment_filter_script = config['DEPENDENCIES']['path_to_alignment_filter_script']
# programs_for_HGT_prediction = config['DEPENDENCIES']['programs_for_HGT_prediction'].split(' ')

############################################### Define folder/file name ################################################

# import modules for tree plot
if plot_tree == 1:
    from ete3 import TreeStyle, NodeStyle, TextFace
    from PIL import Image, ImageDraw, ImageFont

wd = os.getcwd()
op_folder = 'output_ip%s_al%sbp_c%s_e%sbp' % (str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(ending_match_length))

tree_folder =                'tree_folder'
tree_plots_folder =          'tree_plots'
gene_tree_newick_folder =    'gene_tree_newick'
species_tree_folder_plot =   'species_tree_plot'
species_tree_folder_ranger = 'species_tree'
ranger_inputs_folder_name =  'Ranger_input'
ranger_outputs_folder_name = 'Ranger_output'
candidates_file_name =       'HGT_candidates.txt'
candidates_seq_file_name =   'HGT_candidates.fasta'

candidates_file_name_ET =    'HGT_candidates_ET.txt'
candidates_file_name_ET_validated = 'HGT_candidates_ET_validated.txt'
candidates_file_name_ET_validated_fasta = 'HGT_candidates_ET_validated.fasta'

ranger_wd_name =             'Ranger-DTL_wd'
output_tree_folder_name =    'Explicit_tree_output'
tree_image_folder_name =     'combined_tree_images'
ffn_file =                   'combined.ffn'
pwd_prokka =                     '%s/%s'         % (wd, prokka_output)
pwd_grouping_file =              '%s/%s'         % (wd, grouping_file)
pwd_ffn_file =                   '%s/%s'         % (wd, ffn_file)
pwd_ortholog_group_folder =      '%s/%s'         % (wd, ortholog_group_folder_name)
pwd_candidates_file =            '%s/%s/%s'      % (wd, op_folder, candidates_file_name)
pwd_candidates_seq_file =        '%s/%s/%s'      % (wd, op_folder, candidates_seq_file_name)
pwd_candidates_file_ET =         '%s/%s/%s'      % (wd, op_folder, candidates_file_name_ET)
pwd_candidates_file_ET_validated ='%s/%s/%s'      % (wd, op_folder, candidates_file_name_ET_validated)
pwd_candidates_file_ET_validated_fasta ='%s/%s/%s'      % (wd, op_folder, candidates_file_name_ET_validated_fasta)
pwd_ranger_wd =                  '%s/%s/%s/%s/'  % (wd, op_folder, output_tree_folder_name, ranger_wd_name)
pwd_ranger_inputs_folder =       '%s/%s'         % (pwd_ranger_wd, ranger_inputs_folder_name)
pwd_ranger_outputs_folder =      '%s/%s'         % (pwd_ranger_wd, ranger_outputs_folder_name)
pwd_op_tree_folder =             '%s/%s/%s'      % (wd, op_folder, output_tree_folder_name)
pwd_tree_image_folder =          '%s/%s'         % (pwd_op_tree_folder, tree_image_folder_name)
pwd_tree_folder =                '%s/%s'         % (pwd_op_tree_folder, tree_folder)
pwd_tree_plots_folder =          '%s/%s'         % (pwd_op_tree_folder, tree_plots_folder)
pwd_gene_tree_newick_folder =    '%s/%s'         % (pwd_op_tree_folder, gene_tree_newick_folder)
pwd_species_tree_folder_plot =   '%s/%s'         % (pwd_op_tree_folder, species_tree_folder_plot)
pwd_species_tree_folder_ranger = '%s/%s'         % (pwd_op_tree_folder, species_tree_folder_ranger)

# check whether file exist
unfound_inputs = []
for each_input in [pwd_prokka, pwd_grouping_file, pwd_ortholog_group_folder, pwd_phylo_hmm]:
    if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
        unfound_inputs.append(each_input)
if len(unfound_inputs) > 0:
    for each_unfound in unfound_inputs:
        print('%s not found' % each_unfound)
    exit()

########################################################################################################################

# get list of match pair list
candidates_file = open(pwd_candidates_file)
candidates_list = []
for match_group in candidates_file:
    if not match_group.startswith('Recipient'):
        match_group_split = match_group.strip().split('\t')[:2]
        candidates_list.append(match_group_split)


# get all ortholog groups
clusters_original = [os.path.basename(file_name) for file_name in glob.glob('%s/*.fna' % pwd_ortholog_group_folder)]
# remove " ' " from ortholog group name
clusters = []
for cluster_o in clusters_original:
    if "\'" in cluster_o:
        cmd_line_name_fna = cluster_o.replace('\'', '\\\'')
        #cmd_line_name_faa = cmd_line_name_fna.replace('.fna', '.faa')
        cluster_new_fna = cluster_o.replace('\'', '')
        #cluster_new_faa = cluster_new_fna.replace('.fna', '.faa')
        clusters.append(cluster_new_fna)
        os.system('mv %s/%s %s/%s' % (pwd_ortholog_group_folder, cmd_line_name_fna, pwd_ortholog_group_folder, cluster_new_fna))
        #os.system('mv %s/%s %s/%s' % (pwd_ortholog_group_folder, cmd_line_name_faa, pwd_ortholog_group_folder, cluster_new_faa))
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
grouping_profile = open(pwd_grouping_file)
bin_record_list = []
genome_name_list = []
name_to_group_number_dict = {}
name_to_group_number_without_underscore_dict = {}
bin_group_list = []
bin_group_without_underscore_list = []
for each_bin in grouping_profile:
    each_bin_split = each_bin.strip().split(',')
    bin_group = each_bin_split[0]
    bin_group_without_underscore = bin_group.split('_')[0] + bin_group.split('_')[1]
    bin_name = each_bin_split[1]
    name_to_group_number_dict[bin_name] = bin_group
    name_to_group_number_without_underscore_dict[bin_name] = bin_group_without_underscore
    bin_record = BinRecord(bin_name, bin_group, bin_group_without_underscore)
    bin_record_list.append(bin_record)
    genome_name_list.append(bin_name)
    bin_group_list.append(bin_group)
    bin_group_without_underscore_list.append(bin_group_without_underscore)


################################################### Create folders #####################################################

# # create output_tree_folder, gene_tree_seq and gene_tree_txt folder
# if not os.path.isdir(pwd_op_tree_folder):
#     os.mkdir(pwd_op_tree_folder)
#     os.mkdir(pwd_tree_folder)
#     if plot_tree == 1:
#         os.mkdir(pwd_tree_plots_folder)
# else:
#     shutil.rmtree(pwd_op_tree_folder)
#     shutil.rmtree(pwd_op_tree_folder)
#     os.mkdir(pwd_op_tree_folder)
#     os.mkdir(pwd_tree_folder)
#     if plot_tree == 1:
#         os.mkdir(pwd_tree_plots_folder)
#
# # prepare Ranger-DTL working directory
# if not os.path.exists(pwd_ranger_wd):
#     os.makedirs(pwd_ranger_wd)
#     os.makedirs(pwd_ranger_inputs_folder)
#     os.makedirs(pwd_ranger_outputs_folder)
# else:
#     shutil.rmtree(pwd_ranger_wd)
#     shutil.rmtree(pwd_ranger_wd)
#     os.makedirs(pwd_ranger_wd)
#     os.makedirs(pwd_ranger_inputs_folder)
#     os.makedirs(pwd_ranger_outputs_folder)


####################################################### Main Code ######################################################

# get species tree for all input genomes
print('Building species tree for all input genomes')
get_species_tree_wd = '%s/%s' % (wd, 'get_species_tree_wd')

if (os.path.isdir(get_species_tree_wd)) and (len(os.listdir(get_species_tree_wd)) > 1):
    print('get_species_tree_wd folder detected, will skip the alignment step')
else:
    #os.system('rm -r ' + get_species_tree_wd)
    #os.mkdir(get_species_tree_wd)
    get_species_tree_alignment(get_species_tree_wd, pwd_prokka, pwd_phylo_hmm, pwd_hmmsearch_exe, pwd_mafft_exe)
    all_input_genomes = [i for i in os.listdir(pwd_prokka) if os.path.isdir(pwd_prokka + '/' + i)]
    get_species_tree_newick(get_species_tree_wd, all_input_genomes, pwd_fasttree_exe, wd, 'species_tree_all')

candidate_2_predictions_dict = {}
candidate_2_possible_direction_dict = {}

# get gene tree for each orthologous and run Ranger-DTL
print('Get gene tree for each orthologous and run Ranger-DTL')
n = 1
for each_candidates in candidates_list:
    process_name = '___'.join(each_candidates)
    print("Processing %sth of %s candidates: %s" % (str(n), str(len(candidates_list)), process_name))

    # get ortholog_list for each match pairs
    ortholog_list = [] # ortholog_list == gene_member, need to modify!!!!!
    for each in clusters_dict:
        if (each_candidates[0] in clusters_dict[each]) or (each_candidates[1] in clusters_dict[each]):
            ortholog_list += clusters_dict[each]
    if len(ortholog_list) == 0:
        print('No orthologes for the current HGT candidates')
        candidate_2_predictions_dict[process_name] = []
        candidate_2_possible_direction_dict[process_name] = []
    else:
        # get bin name list from candidates gene list
        each_candidates_bin_name = []
        for each_g in each_candidates:
            each_g_new = '_'.join(each_g.split('_')[:-1])
            each_candidates_bin_name.append(each_g_new)

        # uniq ortholog_list, why??????
        gene_member = []
        genome_subset = []
        for each_g in ortholog_list:
            each_g_genome = '_'.join(each_g.split('_')[:-1])
            if each_g not in gene_member:
                gene_member.append(each_g)
            if each_g_genome not in genome_subset:
                genome_subset.append(each_g_genome)
        # get sequences of othorlog group to build gene tree
        seq_file_name_prefix = '___'.join(each_candidates) + '_gene_tree'
        seq_file_name = '%s.seq' % seq_file_name_prefix
        pwd_seq_file = '%s/%s' % (pwd_tree_folder, seq_file_name)
        # output_handle = open(pwd_seq_file, "w")
        # for seq_record in SeqIO.parse(pwd_ffn_file, 'fasta'):
        #     if seq_record.id in gene_member:
        #         SeqIO.write(seq_record, output_handle, 'fasta')
        # output_handle.close()


############################################ get gene tree and species tree ############################################

        # run mafft
        pwd_seq_file_1st_aln = '%s/%s.aln' % (pwd_tree_folder, seq_file_name_prefix)
        #os.system('%s --quiet --maxiterate 1000 --globalpair %s > %s'% (pwd_mafft_exe, pwd_seq_file, pwd_seq_file_1st_aln))

        # run fasttree
        pwd_gene_tree_newick = '%s/%s.newick' % (pwd_tree_folder, seq_file_name_prefix)
        #os.system('%s -nt -quiet %s > %s' % (pwd_fasttree_exe, pwd_seq_file_1st_aln, pwd_gene_tree_newick))

        # plot gene tree
        if plot_tree == 1:
            plot_gene_tree(pwd_gene_tree_newick, 'Gene Tree', process_name, process_name + '_gene_tree', each_candidates, pwd_tree_folder)

        # get species tree subset
        species_tree_file_name = '___'.join(each_candidates) + '_species_tree'
        #get_species_tree_newick(get_species_tree_wd, genome_subset, pwd_fasttree_exe, pwd_tree_folder, species_tree_file_name)

        # plot species tree
        if plot_tree == 1:
            pwd_species_tree_newick_file = '%s/%s.newick' % (pwd_tree_folder, species_tree_file_name)
            plot_species_tree(pwd_species_tree_newick_file, 'Species Tree', process_name, process_name + '_species_tree', each_candidates_bin_name, pwd_tree_folder)


##################################################### Run Ranger-DTL ###################################################

        # define Ranger-DTL input file name
        ranger_inputs_file_name = process_name + '.txt'
        ranger_outputs_file_name = process_name + '_ranger_output.txt'
        pwd_ranger_inputs = '%s/%s' % (pwd_ranger_inputs_folder, ranger_inputs_file_name)
        pwd_ranger_outputs = '%s/%s' % (pwd_ranger_outputs_folder, ranger_outputs_file_name)
        #ranger_inputs_file = open(pwd_ranger_inputs, 'w')

        # read in species tree
        pwd_species_tree_newick = '%s/%s_species_tree.newick' % (pwd_tree_folder, process_name)
        #species_tree = Tree(pwd_species_tree_newick, format=1)
        #species_tree.resolve_polytomy(recursive=True) # solving multifurcations

        # read in gene tree
        #pwd_gene_tree_newick = '%s/%s_gene_tree.newick' % (pwd_tree_folder, process_name)
        #gene_tree = Tree(pwd_gene_tree_newick, format=1)
        #gene_tree.resolve_polytomy(recursive=True) # solving multifurcations

        # change gene tree leaf name for Ranger-DTL
        # for each_gr_leaf in gene_tree:
        #     each_gr_leaf_name_split = each_gr_leaf.name.split('_')
        #     each_gr_leaf.name = '_'.join(each_gr_leaf.name.split('_')[:-1])

        # write species tree and gene tree to Ranger-DTL input file
        #ranger_inputs_file.write('%s\n[&U]%s\n' % (species_tree.write(format=9), gene_tree.write(format=9)))
        #ranger_inputs_file.close()

        # run Ranger-DTL
        ranger_parameters = '-q -D 2 -T 3 -L 1'
        ranger_cmd = '%s %s -i %s -o %s' % (pwd_ranger_exe, ranger_parameters, pwd_ranger_inputs, pwd_ranger_outputs)
        #os.system(ranger_cmd)

        # parse prediction result
        ranger_result = open(pwd_ranger_outputs)
        predicted_transfers = []
        for each_line in ranger_result:
            if 'Transfer' in each_line:
                if not each_line.startswith('The minimum reconciliation cost'):
                    mapping = each_line.strip().split(':')[1].split(',')[1]
                    recipient = each_line.strip().split(':')[1].split(',')[2]
                    donor_p = mapping.split('-->')[1][1:]
                    recipient_p = recipient.split('-->')[1][1:]
                    predicted_transfer = donor_p + '-->' + recipient_p
                    predicted_transfers.append(predicted_transfer)
        candidate_2_predictions_dict[process_name] = predicted_transfers
        # get two possible transfer situation
        candidate_split_group = []
        candidate_split_gene = process_name.split('___')
        candidate_split_gene_only_genome = []
        for each_candidate in candidate_split_gene:
            each_candidate_genome = '_'.join(each_candidate.split('_')[:-1])
            candidate_split_gene_only_genome.append(each_candidate_genome)

        possible_hgt_1 = '%s-->%s' % (candidate_split_gene_only_genome[0], candidate_split_gene_only_genome[1])
        possible_hgt_2 = '%s-->%s' % (candidate_split_gene_only_genome[1], candidate_split_gene_only_genome[0])
        possible_hgts = [possible_hgt_1, possible_hgt_2]
        candidate_2_possible_direction_dict[process_name] = possible_hgts


##################### Combine Species and Gene Tree Together and add Ranger-DTL prediction results #####################

        # read in species/gene tree image
        if plot_tree == 1:
            pwd_species_tree_png = '%s/%s_species_tree.png' % (pwd_tree_folder, process_name)
            pwd_gene_tree_png = '%s/%s_gene_tree.png' % (pwd_tree_folder, process_name)
            species_tree_image = Image.open(pwd_species_tree_png)
            gene_tree_image = Image.open(pwd_gene_tree_png)
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
            #font_arial = ImageFont.truetype("arial.ttf", 24)

            # predicted_transfers
            x_start = 30
            y_start_t = 30 + max(heights)
            y_start = 30 + max(heights) + 40
            #add_text.text((x_start, y_start_t), 'Ranger-DTL predicted HGTs (%s): ' % len(predicted_transfers), (0, 0, 0), font = font_arial)
            add_text.text((x_start, y_start_t), 'Ranger-DTL predicted HGTs (%s): ' % len(predicted_transfers), (0, 0, 0))

            for hgt in predicted_transfers:
                hgt_split = hgt.split('-->')
                hgt_d = hgt_split[0]
                hgt_r = hgt_split[1]

                if x_start < 1800:
                    if hgt in possible_hgts:
                        #add_text.text((x_start, y_start), hgt, (225, 0, 0), font = font_arial)
                        add_text.text((x_start, y_start), hgt, (225, 0, 0))

                        x_start += 200
                    else:
                        if (hgt_d in bin_group_without_underscore_list) and (hgt_r in bin_group_without_underscore_list) and (hgt_d[0] != hgt_r[0]):
                            #add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                            add_text.text((x_start, y_start), hgt, (0, 0, 0))
                            x_start += 200
                        else:
                            #add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                            add_text.text((x_start, y_start), hgt, (0, 0, 0))
                            x_start += 200
                elif x_start >= 1800:
                    x_start = 30
                    y_start += 40
                    if hgt in possible_hgts:
                        #add_text.text((x_start, y_start), hgt, (225, 0, 0), font = font_arial)
                        add_text.text((x_start, y_start), hgt, (225, 0, 0))
                        x_start += 200
                    else:
                        if (hgt_d in bin_group_without_underscore_list) and (hgt_r in bin_group_without_underscore_list) and (hgt_d[0] != hgt_r[0]):
                            #add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                            add_text.text((x_start, y_start), hgt, (0, 0, 0))
                            x_start += 200
                        else:
                            #add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font_arial)
                            add_text.text((x_start, y_start), hgt, (0, 0, 0))
                            x_start += 200

            new_image.save('%s/%s_combined_trees.png' % (pwd_tree_plots_folder, process_name))
            os.remove(pwd_species_tree_png)  # remove species tree
            os.remove(pwd_gene_tree_png)  # remove gene tree
    n += 1

# add results to output file of best blast match approach
print('Add Ranger-DTL predicted direction to HGT_candidates.txt')
combined_output_handle = open(pwd_candidates_file_ET, 'w')
combined_output_validated_handle = open(pwd_candidates_file_ET_validated, 'w')
combined_output_validated_handle.write('Gene_1\tGene_2\tGenome_1_ID\tGenome_2_ID\tIdentity\tEnd_break\tDirection\n' % ())
combined_output_handle.write('Gene_1\tGene_2\tGenome_1_ID\tGenome_2_ID\tIdentity\tEnd_break\tDirection\n' % ())

validated_candidate_list = []
for match_group in open(pwd_candidates_file):
    if not match_group.startswith('Recipient'):
        match_group_split = match_group.strip().split('\t')
        recipient_gene = match_group_split[0]
        donor_gene = match_group_split[1]
        recipient_genome_id = match_group_split[2]
        donor_genome_id = match_group_split[3]
        identity = match_group_split[4]
        end_break = match_group_split[5]
        concatenated = '%s___%s' % (recipient_gene, donor_gene)
        possible_direction = candidate_2_possible_direction_dict[concatenated]
        validated_prediction = 'N/A'
        for each_prediction in candidate_2_predictions_dict[concatenated]:
            if each_prediction in possible_direction:
                validated_prediction = each_prediction

        if (end_break == 'no') and (validated_prediction != 'N/A'):

            if recipient_gene not in validated_candidate_list:
                validated_candidate_list.append(recipient_gene)
            if donor_gene not in validated_candidate_list:
                validated_candidate_list.append(donor_gene)

            combined_output_validated_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_id, donor_genome_id, identity, end_break, validated_prediction))

        combined_output_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (recipient_gene, donor_gene, recipient_genome_id, donor_genome_id, identity, end_break, validated_prediction))

combined_output_handle.close()
combined_output_validated_handle.close()

# export the sequence of validated candidates
combined_output_validated_fasta_handle = open(pwd_candidates_file_ET_validated_fasta, 'w')
for each_candidate in SeqIO.parse(pwd_candidates_seq_file, 'fasta'):
    if each_candidate.id in validated_candidate_list:
        SeqIO.write(each_candidate, combined_output_validated_fasta_handle, 'fasta')
combined_output_validated_fasta_handle.close()


print('\nAll done!')
