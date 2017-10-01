#!/usr/bin/python3

import dendropy
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

# NEEDS, hmmer3, mafft and fasttree

# Call with:    1) Path to folder containing all prokka output folders
#               2) Path to its temporary output folder
#               3) Path to the phylo.hmm file
#               4) Path to write the final alignment and tree
from matplotlib.mlab import phase_spectrum


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def treeMaker (path_to_prokka, path_to_tmp, path_to_hmm, path_to_out):

    # Tests for presence of the tmp folder and deletes it
    if os.path.exists(path_to_tmp):
        os.system('rm -r ' + path_to_tmp)
    os.mkdir(path_to_tmp)

    # List all prokka dirs in the target folder
    prokka_files = [ i for i in os.listdir(path_to_prokka) if os.path.isdir(path_to_prokka + '/' + i) ]

    # Running hmmsearch on each file
    for f in prokka_files:
        # call hmmsearch
        os.system('hmmsearch -o /dev/null --domtblout ' + path_to_tmp + '/' + f + '_hmmout.tbl ' + path_to_hmm + ' ' + path_to_prokka + '/' + f + '/' + f + '.faa')

        # Reading the protein file in a dictionary
        proteinSequence = {}
        for seq_record in SeqIO.parse('%s/%s/%s.faa' % (path_to_prokka, f, f), 'fasta'):
            proteinSequence[seq_record.id] = str(seq_record.seq)

        # Reading the hmmersearch table/extracting the protein part found beu hmmsearch out of the protein/Writing each protein sequence that was extracted to a fasta file (one for each hmm in phylo.hmm
        proteinName = ''
        for line in open('%s/%s_hmmout.tbl' % (path_to_tmp, f), 'r'):
            if line[0] != "#":
                line = re.sub('\s+', ' ', line) # ro remove redundant spaces
                splitLine = line.split(' ')
                if proteinName != splitLine[3]:
                    proteinName = splitLine[3]
                    file_out = open(path_to_tmp + '/' + proteinName + '.fasta', 'a+')
                    x1 = int(splitLine[17]) - 1
                    x2 = int(splitLine[18])
                    export_dna_record(proteinSequence[splitLine[0]][x1: x2], f, '', file_out)
                    file_out.close()

    # Call mafft to align all single fasta files with hmms
    files = os.listdir(path_to_tmp)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    for f in fastaFiles:
        fastaFile1 = path_to_tmp + '/' + f
        fastaFile2 = fastaFile1.replace('.fasta', '_aligned.fasta')
        os.system('mafft --quiet --maxiterate 1000 --globalpair ' + fastaFile1 + ' > ' + fastaFile2 + ' ; rm ' + fastaFile1)

    # concatenating the single alignments
    # create the dictionary
    concatAlignment = {}
    for element in prokka_files:
        concatAlignment[element] = ''
    print(concatAlignment)

    # Reading all single alignment files and append them to the concatenated alignment
    files = os.listdir(path_to_tmp)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    for f in fastaFiles:
        fastaFile = path_to_tmp + '/' + f
        proteinSequence = {}
        print(fastaFile)
        for seq_record_2 in SeqIO.parse(fastaFile, 'fasta'):
            print(seq_record_2)

            proteinName =seq_record_2.id
            proteinSequence[proteinName] = str(seq_record_2.seq)






        alignmentLength = len(proteinSequence[proteinName])

        for element in prokka_files:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    file_out = open(path_to_out + '/alignment.fasta', 'w')
    for element in prokka_files:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')

    # calling fasttree for tree calculation
    os.system('fasttree ' + path_to_out + '/alignment.fasta > ' + path_to_out + '/phylogenticTree.phy')

    # Decomment the two following lines if tree is rooted but should be unrooted
    #phyloTree = dendropy.Tree.get(path='phylogenticTree.phy', schema='newick', rooting='force-unrooted')
    #dendropy.Tree.write_to_path(phyloTree, 'phylogenticTree_unrooted.phy', 'newick')

    #os.system('rm -r ' + path_to_tmp)

treeMaker('prokka', 'hmm_tmp', 'phylo.hmm', '.')
