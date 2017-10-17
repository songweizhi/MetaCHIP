from ete3 import Tree

species_tree = Tree('/Users/songweizhi/Desktop/species_tree_3_rpsO.txt', format = 1)
for each_g_leaf in species_tree:
    print(each_g_leaf)


#     # import species tree and gene tree
#     ranger_species_tree = Tree('%s/species_tree_ranger_%s.txt' % (pwd_species_tree_folder_ranger, candidate), format = 1)
#     ranger_gene_tree = Tree('%s/gene_tree_%s.txt' % (pwd_gene_tree_folder_ranger, candidate), format=1)
#     species_tree_for_plot = '%s/species_tree_plot_%s.txt' % (pwd_species_tree_folder_plot, candidate)
#     # assign new leaf name to gene tree
#     for each_g_leaf in ranger_gene_tree:
#         each_g_leaf_name_split = each_g_leaf.name.split('_')
#         group_number = each_g_leaf_name_split[0]
#         gene_id = each_g_leaf_name_split[1]
#         for each_bin_record in bin_record_list:
#             if group_number == each_bin_record.group_without_underscore:
#                 new_g_leaf_name = '%s_%s_%s' % (group_number, each_bin_record.name, gene_id)
#                 each_g_leaf.name = new_g_leaf_name