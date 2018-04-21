import os
from datetime import datetime
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from ete3 import Tree, ClusterTree



wd = '/Users/songweizhi/Desktop/get_distance_distribution'
tree_plotter_R = '~/R_scripts/newick_tree/tree_plotter.R'
os.chdir(wd)

aln_file = 'species_tree3.aln'
tree_file_out = 'species_tree3.newick'


# read in tree file
aln_handle = AlignIO.read(aln_file, "fasta")

# get distance matrix
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get distance matrix')
calculator = DistanceCalculator('identity')
# calculator = DistanceCalculator('blosum62')
distance_matrix = calculator.get_distance(aln_handle)
print(distance_matrix)
# seq_id_1 = 'bin159'
# seq_id_2 = 'bin1097'
# seq_name_list = distance_matrix.names
# distance_between_seqs = distance_matrix[seq_id_1, seq_id_2]

# construct a tree
# see the difference between UPGMA, NJ and Maximum likelihood
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Construct tree')
constructor = DistanceTreeConstructor()
tree_NJ = constructor.nj(distance_matrix)

# turn to newick format and write out
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Turn to newick format and write out')
tree_NJ_newick = Tree(tree_NJ.format('newick'),format=3)
tree_NJ_newick.write(format=3, outfile=tree_file_out)

# plot the tree
os.system('Rscript %s -t %s' % (tree_plotter_R, tree_file_out))

# cluster distance matrix
