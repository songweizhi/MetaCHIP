import os

wd = '/Users/songweizhi/Desktop/o34'
protein_id_cog_file = 'HGT_candidates_PG_aa_COG_results/protein-id_cog.txt'
cog_stats_file = 'HGT_candidates_PG_aa_COG_results/cog_stats.txt'

pwd_normal_plot_folder = 'Flanking_region_plots/1_Plots_normal_PG_validated'
pwd_normal_plot_folder_high_iden = 'Flanking_region_plots/1_Plots_normal_PG_validated/high_identity'


os.chdir(wd)


# read COG annotatio into dict
Gene_COG_dict = {}
for each_gene_cog in open(protein_id_cog_file):
    each_cog_split = each_gene_cog.strip().split()
    Gene_COG_dict[each_cog_split[0]] = each_cog_split[1]


# read COG annotatio into dict
COG_id_fun_dict = {}
for each_COG in open(cog_stats_file):
    each_COG_split = each_COG.strip().split('\t')
    each_COG_id = each_COG_split[0]
    each_COG_fun = each_COG_split[1]
    COG_id_fun_dict[each_COG_id] = each_COG_fun


# get gene list
HGT_gene_list = []
HGT_gene_pair_list = []
for each in open('HGT_candidates_PG_validated.txt'):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        Gene_1 = each_split[0]
        Gene_2 = each_split[1]
        Identity = float(each_split[4])
        End_break = each_split[5]
        Direction = each_split[7]
        if Identity >= 95:
            if Gene_1 not in HGT_gene_list:
                HGT_gene_list.append(Gene_1)
            if Gene_2 not in HGT_gene_list:
                HGT_gene_list.append(Gene_2)
            HGT_gene_pair_list.append([Gene_1, Gene_2])


# get COG ID and function for gene list
for each_gene in sorted(HGT_gene_list):

    # get COG ID
    gene_cog = 'NA'
    if each_gene in Gene_COG_dict:
        gene_cog = Gene_COG_dict[each_gene]

    # get COG function
    gene_fun = 'NA'
    if gene_cog in COG_id_fun_dict:
        gene_fun = COG_id_fun_dict[gene_cog]

    for_print = '%s\t%s\t%s' % (each_gene, gene_cog, gene_fun)
    print(for_print)


for each_gene_pair in HGT_gene_pair_list:
    print(each_gene_pair)
    possible_plot_1 = '%s___%s.eps' % (each_gene_pair[0], each_gene_pair[1])
    possible_plot_2 = '%s___%s.eps' % (each_gene_pair[1], each_gene_pair[0])
    pwd_possible_plot_1 = '%s/%s' % (pwd_normal_plot_folder, possible_plot_1)
    pwd_possible_plot_2 = '%s/%s' % (pwd_normal_plot_folder, possible_plot_2)

    if os.path.isfile(pwd_possible_plot_1) == True:
        os.system('cp %s %s/' % (pwd_possible_plot_1, pwd_normal_plot_folder_high_iden))

    if os.path.isfile(pwd_possible_plot_2) == True:
        os.system('cp %s %s/' % (pwd_possible_plot_2, pwd_normal_plot_folder_high_iden))





