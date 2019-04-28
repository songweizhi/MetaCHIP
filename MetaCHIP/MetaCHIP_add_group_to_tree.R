# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

# Usage example
packages<-c("ape", "tools", "optparse")
invisible(suppressMessages(check.packages(packages)))

# install.packages('ape')
# suppressWarnings(library(ape))
# suppressWarnings(library(tools))
# suppressWarnings(library(optparse))

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
              metavar="character"),
  
  make_option(c("-f", "--label_font"), 
              type="double", 
              help="label fontsize on the tree plot", 
              metavar="double"),

  make_option(c("-s", "--label_shift"),
              type="double",
              help="label shift on the tree plot",
              metavar="double"),
  
  make_option(c("-v", "--tree_height"),
              type="double",
              help="height of tree plot",
              metavar="double"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
wd = getwd()
grouping_file = opt$grouping
tree_file_in = opt$tree
label_shift = opt$label_shift
label_font = opt$label_font
tree_height = opt$tree_height

tree_file_in_name_no_extension = file_path_sans_ext(basename(grouping_file))
tree_txt_file_with_group = paste(tree_file_in_name_no_extension, 'tree.txt', sep = '_')
tree_txt_file_only_group = paste(tree_file_in_name_no_extension, 'tree_group_only.txt', sep = '_')
tree_plot_file_with_group = paste(tree_file_in_name_no_extension, 'tree.jpg', sep = '_')
tree_plot_file_only_group = paste(tree_file_in_name_no_extension, 'tree_group_only.jpg', sep = '_')
tree_plot_file_only_group_u = paste(tree_file_in_name_no_extension, 'tree_group_only_unrooted.jpg', sep = '_')

pwd_tree_txt_file_with_group = paste(wd, tree_txt_file_with_group, sep = '/')
pwd_tree_txt_file_only_group = paste(wd, tree_txt_file_only_group, sep = '/')
pwd_tree_plot_file_with_group = paste(wd, tree_plot_file_with_group, sep = '/')
pwd_tree_plot_file_only_group = paste(wd, tree_plot_file_only_group, sep = '/')
pwd_tree_plot_file_only_group_u = paste(wd, tree_plot_file_only_group_u, sep = '/')

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
jpeg(pwd_tree_plot_file_with_group, width = 2000, height = tree_height, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_with_group, 'phylogram', font = 1, cex = label_font, label.offset = label_shift, align.tip.label = 2, lab4ut = 'axial')
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
jpeg(pwd_tree_plot_file_only_group, width = 2000, height = tree_height, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_only_group, 'phylogram', font = 1, cex = label_font, label.offset = label_shift, align.tip.label = 2, lab4ut = 'axial')
dev.off()
#write.tree(SCG_tree_only_group, file=pwd_tree_txt_file_only_group)


# plot SCG tree only group (unrooted)
SCG_tree_only_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_only_group$tip.label)) {
  label_name = SCG_tree_only_group$tip.label[i]
  label_name_row_num = which(grouping_df$V2 == label_name)
  group_id = grouping_df$V1[label_name_row_num]
  SCG_tree_only_group$tip.label[i] = paste(group_id)
  i = i + 1}
# plot tree only group
jpeg(pwd_tree_plot_file_only_group_u, width = 2000, height = 2000, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_only_group, 'u', font = 1, cex = label_font, label.offset = label_shift, lab4ut = 'axial')
dev.off()

