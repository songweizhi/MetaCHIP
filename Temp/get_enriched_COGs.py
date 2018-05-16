
usage = '''
python /Users/songweizhi/PycharmProjects/MetaCHIP/get_AR_COGs.py
'''

# annotation resultpwds
cog_annotation_results = '/Users/songweizhi/Desktop/protein-id_cog.txt'
enriched_cogs = ['O', 'C', 'I']

# db files
whog_file = '/Users/songweizhi/Softwares/COG_annotation_db/whog'

# get cog_id to description and category dict
cog_id_description_dict = {}
cog_id_category_dict = {}
for each_cog in open(whog_file):
    if each_cog.startswith('['):
        each_cog_split = each_cog.strip().split(' ')
        each_cog_category = each_cog_split[0][1:-1]
        each_cog_id = each_cog_split[1]
        each_cog_description = ' '.join(each_cog_split[2:])
        cog_id_description_dict[each_cog_id] = each_cog_description
        cog_id_category_dict[each_cog_id] = each_cog_category

for each_gene in open(cog_annotation_results):
    each_gene_split = each_gene.strip().split('\t')
    gene_id = each_gene_split[0]
    gene_cog = each_gene_split[1]
    gene_cate = each_gene_split[2]

    for each_enriched_cog in enriched_cogs:
        if each_enriched_cog in gene_cate:
            print('%s\t%s\t%s\t%s' % (gene_cate, gene_cog, gene_id, cog_id_description_dict[gene_cog]))
