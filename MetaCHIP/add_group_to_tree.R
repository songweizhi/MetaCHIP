#!/usr/bin/env Rscript

######################################## Usage ########################################

# usgae
# Rscript add_group_to_tree.R -t input_tree.newick -g grouping_file.txt

#######################################################################################


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
tree_txt_file_with_group = paste(tree_file_in_name_no_extension, 'with_group.newick', sep = '_')
tree_txt_file_only_group = paste(tree_file_in_name_no_extension, 'only_group.newick', sep = '_')

pwd_tree_txt_file_with_group = paste(wd, tree_txt_file_with_group, sep = '/')
pwd_tree_txt_file_only_group = paste(wd, tree_txt_file_only_group, sep = '/')

# read in grouping file
grouping_df = read.csv(grouping_file, header = FALSE)


#################### get tree with group ####################

SCG_tree_with_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_with_group$tip.label)) {
  label_name = SCG_tree_with_group$tip.label[i]
  label_name_row_num = which(grouping_df$V2 == label_name)
  group_id = grouping_df$V1[label_name_row_num]
  SCG_tree_with_group$tip.label[i] = paste(group_id, SCG_tree_with_group$tip.label[i], sep = '_')
  i = i + 1}

# write out tree
write.tree(SCG_tree_with_group, file=pwd_tree_txt_file_with_group)


#################### get tree with group only ####################

SCG_tree_only_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_only_group$tip.label)) {
  label_name = SCG_tree_only_group$tip.label[i]
  label_name_row_num = which(grouping_df$V2 == label_name)
  group_id = grouping_df$V1[label_name_row_num]
  SCG_tree_only_group$tip.label[i] = paste(group_id)
  i = i + 1}

# write out tree
write.tree(SCG_tree_only_group, file=pwd_tree_txt_file_only_group)


