
import os
import glob
from Bio import SeqIO
import matplotlib.pyplot as plt


# usearch_parser
UCLUST_output = '/Users/songweizhi/Desktop/mNC_clusters_0.3.uc.t5'
combined_faa = '/Users/songweizhi/Desktop/mNC_combined.faa'
min_gene_num = 3




# srote clustering results into dict
cluster_to_gene_member_dict = {}
for each in open(UCLUST_output):
    each_split = each.strip().split('\t')
    cluster_id = each_split[1]
    gene_member = each_split[8]
    if cluster_id not in cluster_to_gene_member_dict:
        cluster_to_gene_member_dict[cluster_id] = {gene_member}
    else:
        cluster_to_gene_member_dict[cluster_id].add(gene_member)
#print('UCLUST\tThe total number of predicted clusters: %s' % len(cluster_to_gene_member_dict))


# get qualified clusters
qualified_clusters = {}
gene_list_counted = set()
for cluster in cluster_to_gene_member_dict:
    if len(cluster_to_gene_member_dict[cluster]) >= min_gene_num:
        qualified_clusters[cluster] = cluster_to_gene_member_dict[cluster]

        for gene in cluster_to_gene_member_dict[cluster]:
            gene_list_counted.add(gene)


print('UCLUST\tThe number of clusters with at least %s genes: %s' % (min_gene_num, len(qualified_clusters)))
print('UCLUST\tThe number of genes included in qualified clusters: %s' % (len(gene_list_counted)))
print('UCLUST\tAverage gene number per cluster: %s' % float("{0:.1f}".format((len(gene_list_counted)/len(qualified_clusters)))))
print()


# plot cluster size distribution for UCLUST
gene_number_list = []
for cluster in qualified_clusters:
    gene_number_list.append(len(qualified_clusters[cluster]))


# plot identity distribution
# plt.hist(gene_number_list, bins=60, linewidth=1, alpha=0.5)
# plt.show()


files = '/Users/songweizhi/Desktop/mNC_homologues/*.faa'
file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]
cluster_to_gene_member_dict_get_homo = {}
gene_list_counted_get_homo = set()
for file in file_list:
    pwd_file = '/Users/songweizhi/Desktop/mNC_homologues/%s' % file
    cluster_to_gene_member_dict_get_homo[file] = set()
    for gene in SeqIO.parse(pwd_file, 'fasta'):
        gene_id = str(gene.id)[3:]
        cluster_to_gene_member_dict_get_homo[file].add(gene_id)
        gene_list_counted_get_homo.add(gene_id)


gene_number_list_get_homo = []
for cluster in cluster_to_gene_member_dict_get_homo:
    gene_number_list_get_homo.append(len(cluster_to_gene_member_dict_get_homo[cluster]))


# plot identity distribution
# plt.hist(gene_number_list_get_homo, bins=80, linewidth=1, alpha=0.5)
# plt.show()


print('GET_HO\tThe number of clusters with at least %s genes: %s' % (min_gene_num, len(cluster_to_gene_member_dict_get_homo)))
print('GET_HO\tThe number of genes included in qualified clusters: %s' % (len(gene_list_counted_get_homo)))
print('GET_HO\tAverage gene number per cluster: %s' % float("{0:.1f}".format((len(gene_list_counted_get_homo)/len(cluster_to_gene_member_dict_get_homo)))))




