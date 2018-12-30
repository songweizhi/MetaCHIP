import os

wd = '/Users/songweizhi/Desktop/recentHGT_c'
HGT_PG = 'HGT_candidates_PG_validated.txt'
protein_id_cog_file = 'protein-id_cog.txt'
cog_stats_file = 'cog_stats.txt'


os.chdir(wd)


recent_HGT_dict = {}
for each in open(HGT_PG):
    each_split = each.strip().split('\t')
    if not each.startswith('Gene_1'):
        identity = float(each_split[4])
        if identity >= 99:
            recent_HGT_dict['%s,%s' % (each_split[0], each_split[1])] = identity


Gene_COG_dict = {}
for each_gene_cog in open(protein_id_cog_file):
    each_cog_split = each_gene_cog.strip().split()
    Gene_COG_dict[each_cog_split[0]] = each_cog_split[1]


COG_id_fun_dict = {}
for each_COG in open(cog_stats_file):
    each_COG_split = each_COG.strip().split('\t')
    each_COG_id = each_COG_split[0]
    each_COG_fun = each_COG_split[1]
    COG_id_fun_dict[each_COG_id] = each_COG_fun


for each_hgt in recent_HGT_dict:

    each_hgt_split = each_hgt.split(',')
    gene_1 = each_hgt_split[0]
    gene_2 = each_hgt_split[1]

    gene_1_cog = 'NA'
    if gene_1 in Gene_COG_dict:
        gene_1_cog = Gene_COG_dict[gene_1]

    gene_2_cog = 'NA'
    if gene_2 in Gene_COG_dict:
        gene_2_cog = Gene_COG_dict[gene_2]

    gene_1_fun = 'NA'
    if gene_1_cog in COG_id_fun_dict:
        gene_1_fun = COG_id_fun_dict[gene_1_cog]

    gene_2_fun = 'NA'
    if gene_2_cog in COG_id_fun_dict:
        gene_2_fun = COG_id_fun_dict[gene_2_cog]

    gene_cog_combined = ''
    gene_fun_combined = ''
    if gene_1_cog == gene_2_cog:
        gene_cog_combined = gene_1_cog
        gene_fun_combined = gene_1_fun
    else:
        gene_cog_combined = '%s|%s' % (gene_1_cog, gene_2_cog)
        gene_fun_combined = '%s|%s' % (gene_1_fun, gene_2_fun)

    for_print = '%s\t%s\t%s\t%s\t%s' % (each_hgt.split(',')[0], each_hgt.split(',')[1], recent_HGT_dict[each_hgt], gene_cog_combined, gene_fun_combined)
    print(for_print)


