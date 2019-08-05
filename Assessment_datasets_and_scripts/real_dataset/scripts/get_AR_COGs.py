import os
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt

usage = '''

python /Users/songweizhi/PycharmProjects/MetaCHIP/Test_datasets_and_scripts/real_dataset/scripts/get_AR_COGs.py

'''

# annotation resultpwds
protein_id_cog_folder = '/Users/songweizhi/Desktop/protein_id_cog_folder'


# db files
type2cog_file = '/Users/songweizhi/Softwares/ARDB/ARDBflatFiles/type2cog.tab'
whog_file = '/Users/songweizhi/Softwares/COG_annotation_db/whog'


Ribosomal_protein_related_COGs = {'COG0092', 'COG2139', 'COG0093', 'COG2007', 'COG0099', 'COG0828', 'COG0094', 'COG5051', 'COG2004', 'COG2238', 'COG0081', 'COG2163', 'COG0052', 'COG2097', 'COG0230', 'COG1890', 'COG1727', 'COG0360', 'COG4352', 'COG2174', 'COG0203', 'COG0048', 'COG0255', 'COG1825', 'COG0268', 'COG1997', 'COG0051', 'COG2051', 'COG0359', 'COG0244', 'COG2058', 'COG4830', 'COG0257', 'COG0089', 'COG2053', 'COG0200', 'COG2125', 'COG0100', 'COG0098', 'COG0080', 'COG2451', 'COG0184', 'COG1358', 'COG1841', 'COG0292', 'COG0227', 'COG0261', 'COG0256', 'COG4919', 'COG2126', 'COG0090', 'COG0087', 'COG1471', 'COG0539', 'COG0049', 'COG2167', 'COG0228', 'COG0522', 'COG1632', 'COG2075', 'COG0091', 'COG0333', 'COG0222', 'COG0185', 'COG1552', 'COG4901', 'COG0267', 'COG5045', 'COG1998', 'COG1383', 'COG0335', 'COG0197', 'COG0186', 'COG0291', 'COG0238', 'COG2157', 'COG0096', 'COG1911', 'COG0088', 'COG1631', 'COG0198', 'COG0097', 'COG0102', 'COG0254', 'COG0211', 'COG1717', 'COG0199', 'COG2147', 'COG0103'}
Ribosomal_protein_related_COGs = set()


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



files = '%s/*.txt' % protein_id_cog_folder
file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]


percent_list = []
for each_file in sorted(file_list):

    pwd_each_file = '%s/%s' % (protein_id_cog_folder, each_file)

    num = 0
    total_gene = 0
    for each_gene in open(pwd_each_file):
        gene_id = each_gene.strip().split('\t')[0]
        gene_cog_id = each_gene.strip().split('\t')[1]
        #print(gene_cog_id)
        if gene_cog_id not in Ribosomal_protein_related_COGs:
            total_gene += 1
            if gene_cog_id in ar_related_cogs:
                print('%s\t%s\t%s' % (gene_id, gene_cog_id, cog_id_description_dict[gene_cog_id]))
                num += 1


    #print('Number: %s' % str(num))
    print('Percent: %s\t%s' % (each_file, str(num/total_gene)))

    percent_list.append(num/total_gene)


# # Multiple box plots on one Axes
# fig, ax = plt.subplots()
# ax.boxplot(percent_list)
# plt.show()
