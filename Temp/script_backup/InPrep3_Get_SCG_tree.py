#!/usr/bin/python
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def treeMaker(path_to_prokka, path_to_hmm, pwd_hmmsearch_exe, pwd_mafft_exe, pwd_fasttree_exe, plot_tree):

    # Tests for presence of the tmp folder and deletes it
    tmp_folder = 'get_species_tree_wd'
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
        #os.system('hmmsearch -o /dev/null --domtblout %s/%s_hmmout.tbl %s %s/%s/%s.faa' % (tmp_folder, f, path_to_hmm, path_to_prokka, f, f))
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
    print('Running mafft...')
    for f in fastaFiles:
        fastaFile1 = '%s/%s' % (tmp_folder, f)
        fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
        os.system(pwd_mafft_exe + ' --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)

    # concatenating the single alignments
    # create the dictionary
    print('Concatenating alignments...')
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
            proteinName =seq_record_2.id
            proteinSequence[proteinName] = str(seq_record_2.seq)
            alignmentLength = len(proteinSequence[proteinName])

        for element in prokka_files:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    file_out = open('./species_tree.aln', 'w')
    for element in prokka_files:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()

    # calling fasttree for tree calculation
    print('Running fasttree...')
    os.system('%s -quiet species_tree.aln > species_tree.newick' % pwd_fasttree_exe)

    # Decomment the two following lines if tree is rooted but should be unrooted
    #phyloTree = dendropy.Tree.get(path='phylogenticTree.phy', schema='newick', rooting='force-unrooted')
    #dendropy.Tree.write_to_path(phyloTree, 'phylogenticTree_unrooted.phy', 'newick')

    # plot species tree
    if plot_tree == 1:
        print('Plot species tree')

        tree = Tree('species_tree.newick', format=1)
        # set tree parameters
        ts = TreeStyle()
        ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
        ts.show_leaf_name = 0
        # set tree title text parameters
        ts.title.add_face(TextFace('Species_Tree',
                                   fsize=8,
                                   fgcolor='black',
                                   ftype='Arial',
                                   tight_text=False),
                          column=0)  # tree title text setting
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
                ns["shape"] = "circle"  # dot shape: circle, square or sphere
                ns["size"] = 0  # dot size
                ns['hz_line_width'] = 0.5  # branch line width
                ns['vt_line_width'] = 0.5  # branch line width
                ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
                ns['vt_line_type'] = 0  # branch line type
                ns["fgcolor"] = "blue"  # the dot setting
                each_node.add_face(TextFace(each_node.name,
                                            fsize=5,
                                            fgcolor='black',
                                            tight_text=False,
                                            bold=False),
                                   column=0,
                                   position='branch-right')  # leaf node the node name text setting

                each_node.set_style(ns)

            # non-leaf node parameters
            else:
                nlns = NodeStyle()
                nlns["size"] = 0  # dot size
                # nlns["rotation"] = 45
                each_node.add_face(TextFace(each_node.name,
                                            fsize=3,
                                            fgcolor='black',
                                            tight_text=False,
                                            bold=False),
                                   column=5,
                                   position='branch-top')  # non-leaf node name text setting)

                each_node.set_style(nlns)

        tree.render('species_tree' + '.png', w=900, units="px", tree_style=ts)  # set figures size

    if plot_tree == 0:
        print('The built species tree was exported to species_tree.newick')
    else:
        print('The built species tree was exported to species_tree.newick and species_tree.png')


parser = argparse.ArgumentParser()

parser.add_argument('-prokka_output',
                    required=True,
                    help='the folder holds Prokka outputs for input genomes')

parser.add_argument('-hmm',
                    required=True,
                    help='the phylo.hmm file')

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

parser.add_argument('-plot_tree',
                    required=False,
                    action='store_true',
                    help='to plot species tree')

args = vars(parser.parse_args())
prokka_output = args['prokka_output']
hmm_file = args['hmm']
pwd_hmmsearch_exe = args['hmmsearch']
pwd_mafft_exe = args['mafft']
pwd_fasttree_exe = args['fasttree']
plot_tree = args['plot_tree']


if plot_tree == 1:
    from ete3 import Tree, TreeStyle, NodeStyle, TextFace

treeMaker(prokka_output, hmm_file, pwd_hmmsearch_exe, pwd_mafft_exe, pwd_fasttree_exe, plot_tree)
