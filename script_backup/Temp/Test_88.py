
import os
from datetime import datetime
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from ete3 import Tree, ClusterTree
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster




wd = '/Users/songweizhi/Desktop/get_distance_distribution'
tree_plotter_R = '~/R_scripts/newick_tree/tree_plotter.R'
os.chdir(wd)

aln_file = 'species_tree2.aln'
tree_file_out = 'species_tree2.newick'
taxon_assignment = 'taxon_assignment.txt'


# read in bin id to taxon dict
bin_id_to_taxon_dict = {}
for each in open(taxon_assignment):
    each_split = each.strip().split(',')
    bin_id_to_taxon_dict[each_split[0]] = each_split[1]


# read in tree file
aln_handle = AlignIO.read(aln_file, "fasta")

# get distance matrix
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get distance matrix')
calculator = DistanceCalculator('identity')
# calculator = DistanceCalculator('blosum62')
distance_matrix = calculator.get_distance(aln_handle)
#print(distance_matrix)

# construct a tree
# see the difference between UPGMA, NJ and Maximum likelihood
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Construct tree')
constructor = DistanceTreeConstructor()
tree_NJ = constructor.nj(distance_matrix)

# turn to newick format and write out
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Turn to newick format and write out')
tree_NJ_newick = Tree(tree_NJ.format('newick'),format=3)
for each_leaf in tree_NJ_newick:
    #print(each_leaf.name)
    each_leaf.name = bin_id_to_taxon_dict[each_leaf.name]


tree_NJ_newick.write(format=3, outfile=tree_file_out)

# plot the tree
os.system('Rscript %s -t %s' % (tree_plotter_R, tree_file_out))


# get full matrix
seq_name_list = distance_matrix.names

all_distances_lol = []
for each_seq_1 in seq_name_list:
    current_seq_list = []
    for each_seq_2 in seq_name_list:
        distance_between_seqs = distance_matrix[each_seq_1, each_seq_2]
        distance_between_seqs = float("{0:.1f}".format(distance_between_seqs))
        current_seq_list.append(distance_between_seqs)
    all_distances_lol.append(current_seq_list)

all_distances_lol_array = np.array(all_distances_lol)

print(all_distances_lol_array)


cluster = linkage(all_distances_lol_array, method='weighted')

flat_clusters = fcluster(cluster, 0.75, criterion='distance')

print(flat_clusters)

print(seq_name_list)




out = open('clustering.txt', 'w')
n = 0
for each in seq_name_list:
    #print('%s\t%s\t%s' % (flat_clusters[n], each, bin_id_to_taxon_dict[each]))
    out.write('%s\t%s\t%s\n' % (bin_id_to_taxon_dict[each], flat_clusters[n], each))
    n += 1

out.close()


os.system('cat clustering.txt | sort > clustering_sorted.txt')
# distance_between_seqs = distance_matrix[seq_id_1, seq_id_2]

for each in open('clustering_sorted.txt'):
    print(each.strip())

print('The number of sequences: %s' % len(seq_name_list))
