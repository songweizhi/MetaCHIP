
# annotation results
cog_annotation_results = '/Users/songweizhi/Desktop/COG_and_KO/CF_protein-id_cog.txt'
#cog_annotation_results = '/Users/songweizhi/Desktop/COG_and_KO/NorthSea_protein-id_cog.txt'
#cog_annotation_results = '/Users/songweizhi/Desktop/COG_and_KO/Hospital_protein-id_cog.txt'
#cog_annotation_results = '/Users/songweizhi/Desktop/COG_and_KO/MetaBAT_protein-id_cog.txt'


# db files
type2cog_file = '/Users/songweizhi/Desktop/ARDBflatFiles/type2cog.tab'
whog_file = '/Users/songweizhi/Desktop/whog'

# get AR related COGs
ar_related_cogs = []
for each in open(type2cog_file):
    cog_id = 'COG%s' % each.strip().split('\t')[1]
    if cog_id not in ar_related_cogs:
        ar_related_cogs.append(cog_id)

# get cog_id_description_dict
cog_id_description_dict = {}
for each_cog in open(whog_file):
    if each_cog.startswith('['):
        each_cog_split = each_cog.strip().split(' ')
        each_cog_id = each_cog_split[1]
        each_cog_description = ' '.join(each_cog_split[2:])
        cog_id_description_dict[each_cog_id] = each_cog_description


num = 0
total_gene = 0
for each_gene in open(cog_annotation_results):
    gene_id = each_gene.strip().split('\t')[0]
    gene_cog_id = each_gene.strip().split('\t')[1]
    total_gene += 1
    if gene_cog_id in ar_related_cogs:
        print('%s\t%s\t%s' % (gene_id, gene_cog_id, cog_id_description_dict[gene_cog_id]))
        num += 1

print('Number: %s' % str(num))
print('Percent: %s' % str(num/total_gene))
