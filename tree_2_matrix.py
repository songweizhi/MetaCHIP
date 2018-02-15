import os
import dendropy
from dendropy import treecalc


os.chdir('/Users/songweizhi/Desktop')


tree = dendropy.Tree.get_from_path("species_tree.newick", "newick")
#pdm = treecalc.phylogenetic_distance_matrix(tree)
pdc = tree.phylogenetic_distance_matrix()
print(pdc)
print(type(pdc))



for i, t1 in enumerate(tree.taxon_namespace[:-1]):
    for t2 in tree.taxon_namespace[i+1:]:
        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))





















