from Bio import Phylo


tree_file_in = '/Users/songweizhi/Desktop/bac120_r86.1.tree'
tree_file_out = '/Users/songweizhi/Desktop/bac120_r86.1_updated.tree'

trees = Phylo.read(tree_file_in, 'newick')
print(trees)





#Phylo.write(trees, tree_file_out, 'newick')