suppressWarnings(library(ape))
suppressWarnings(library(tools))
suppressWarnings(library(optparse))

# usgae
# Rscript ~/R_scripts/newick_tree/add_group_to_tree.R -t species_tree.newick -g species_tree_grouping_modified.txt

option_list = list(
  
  make_option(c("-t", "--tree"), 
              type="character", 
              help="tree file", 
              metavar="character"),

  make_option(c("-g", "--grouping"), 
              type="character", 
              help="grouping file", 
              metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
wd = getwd()
grouping_file = opt$grouping
tree_file_in = opt$tree

tree_file_in_name_no_extension = file_path_sans_ext(basename(grouping_file))
tree_txt_file_with_group = paste(tree_file_in_name_no_extension, 'with_group.txt', sep = '_')
tree_txt_file_only_group = paste(tree_file_in_name_no_extension, 'only_group.txt', sep = '_')
tree_plot_file_with_group = paste(tree_file_in_name_no_extension, 'with_group.jpg', sep = '_')
tree_plot_file_only_group = paste(tree_file_in_name_no_extension, 'only_group.jpg', sep = '_')
pwd_tree_txt_file_with_group = paste(wd, tree_txt_file_with_group, sep = '/')
pwd_tree_txt_file_only_group = paste(wd, tree_txt_file_only_group, sep = '/')
pwd_tree_plot_file_with_group = paste(wd, tree_plot_file_with_group, sep = '/')
pwd_tree_plot_file_only_group = paste(wd, tree_plot_file_only_group, sep = '/')

# read in grouping file
grouping_df = read.csv(grouping_file, header = FALSE)

# plot SCG tree with group
SCG_tree_with_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_with_group$tip.label)) {
  label_name = SCG_tree_with_group$tip.label[i]
  label_name_row_num = which(grouping_df$V2 == label_name)
  group_id = grouping_df$V1[label_name_row_num]
  SCG_tree_with_group$tip.label[i] = paste(group_id, SCG_tree_with_group$tip.label[i], sep = '_')
  i = i + 1}
# plot tree with group
jpeg(pwd_tree_plot_file_with_group, width = 1600, height = 1600, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_with_group, 'u', font = 1, cex = 0.7, label.offset = 0.02, lab4ut = 'axial')
dev.off()
#write.tree(SCG_tree_with_group, file=pwd_tree_txt_file_with_group)

# plot SCG tree only group
SCG_tree_only_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_only_group$tip.label)) {
  label_name = SCG_tree_only_group$tip.label[i]
  label_name_row_num = which(grouping_df$V2 == label_name)
  group_id = grouping_df$V1[label_name_row_num]
  SCG_tree_only_group$tip.label[i] = paste(group_id)
  i = i + 1}
# plot tree only group
jpeg(pwd_tree_plot_file_only_group, width = 1600, height = 1600, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_only_group, 'u', font = 1, cex = 0.7, label.offset = 0.02, lab4ut = 'axial')
dev.off()
#write.tree(SCG_tree_only_group, file=pwd_tree_txt_file_only_group)
