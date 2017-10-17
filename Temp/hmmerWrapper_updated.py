#!/usr/bin/python3
import os
import re
#import dendropy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import Phylo

# NEEDS, hmmer3, mafft and fasttree
# Call with:    1) Path to folder containing all prokka output folders
#               2) Path to its temporary output folder
#               3) Path to the phylo.hmm file
#               4) Path to write the final alignment and tree


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def get_species_tree_alignment(tmp_folder, path_to_prokka, pwd_hmmsearch_exe, pwd_mafft_exe):

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
        proteinName = ''
        for line in open('%s/%s_hmmout.tbl' % (tmp_folder, f), 'r'):
            if line[0] != "#":
                line = re.sub('\s+', ' ', line) # ro remove redundant spaces
                splitLine = line.split(' ')
                if proteinName != splitLine[3]:
                    proteinName = splitLine[3]
                    file_out = open(tmp_folder + '/' + proteinName + '.fasta', 'a+')
                    x1 = int(splitLine[17]) - 1
                    x2 = int(splitLine[18])
                    export_dna_record(proteinSequence[splitLine[0]][x1: x2], f, '', file_out)
                    file_out.close()

    # Call mafft to align all single fasta files with hmms
    files = os.listdir(tmp_folder)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    print('Running mafft...')
    for f in fastaFiles:
        fastaFile1 = '%s/%s' % (tmp_folder, f)
        fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
        os.system(pwd_mafft_exe + ' --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)


def get_species_tree_newick(tmp_folder, each_subset, pwd_fasttree_exe):
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
    file_out = open('./alignment.fasta', 'w')
    # for element in prokka_files:
    for element in each_subset:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()

    # calling fasttree for tree calculation
    os.system('%s -quiet alignment.fasta > phylogenticTree.phy' % pwd_fasttree_exe)

    # Decomment the two following lines if tree is rooted but should be unrooted
    # phyloTree = dendropy.Tree.get(path='phylogenticTree.phy', schema='newick', rooting='force-unrooted')
    # dendropy.Tree.write_to_path(phyloTree, 'phylogenticTree_unrooted.phy', 'newick')
    # os.system('rm -r ' + tmp_folder)
    species_tree = Phylo.read('phylogenticTree.phy', 'newick')
    Phylo.draw_ascii(species_tree)


os.chdir('/Users/songweizhi/Desktop/hmmTree')
pwd_hmmsearch_exe = '/Users/songweizhi/Softwares/hmmer/hmmer-3.1b2-macosx-intel/binaries/hmmsearch'
pwd_fasttree_exe = '/Users/songweizhi/Softwares/FastTree/fasttree'
pwd_mafft_exe = 'mafft'
path_to_prokka = 'prokka'
path_to_hmm = 'phylo.hmm'
tmp_folder = 'hmm_tmp'
subset_list = [['AMS', 'AAM', 'AKV', 'AMAU', 'AMAC', 'BDS', 'BAD', 'BGC', 'BNM', 'BHS'],
               ['AMS', 'AAM', 'AKV', 'AMAU', 'AMAC'],
               ['BDS', 'BAD', 'BGC', 'BNM', 'BHS'],
               ['AKV', 'AMAU', 'BDS', 'BAD']]


get_species_tree_alignment(tmp_folder, path_to_prokka, pwd_hmmsearch_exe, pwd_mafft_exe)


for each_subset in subset_list:
    get_species_tree_newick(tmp_folder, each_subset, pwd_fasttree_exe)
