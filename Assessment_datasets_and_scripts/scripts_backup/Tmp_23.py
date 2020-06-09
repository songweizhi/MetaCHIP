from Bio import Phylo


def remove_hyphen_from_branch_length(tree_in, tree_out, tree_format):
    tree_in = Phylo.parse(tree_in, tree_format)
    Phylo.write(tree_in, tree_out, tree_format)


tree_in  = '/Users/songweizhi/Desktop/in/DF_bin4_01620___MilliQ6_bin17_01934.txt'
tree_out = '/Users/songweizhi/Desktop/in/DF_bin4_01620___MilliQ6_bin17_01934._renamed.txt'
remove_hyphen_from_branch_length(tree_in, tree_out, 'newick')



print(float("1e-05"))

print(0.00001)

print(type(float("1e-05")))



print(0.00001 == float("1e-05"))