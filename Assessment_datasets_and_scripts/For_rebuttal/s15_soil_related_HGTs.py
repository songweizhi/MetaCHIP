import os
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def bins_labels(bins, **kwargs):
    bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
    plt.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins, **kwargs)
    plt.xlim(bins[0], bins[-1])


wd = '/Users/songweizhi/Desktop/nature_rebuttal'
os.chdir(wd)


# input files
genome_meadata_file = 'nature10571_genome_metadata_2094.txt'
overall_10255_HGTs = 'Total2094_g664_HGTs_PG_validated.txt'
overall_10255_HGTs_seq = 'Total2094_g664_HGTs_PG_aa.fasta'


# output files
soil_HGT_recent_faa = 'soil_HGT_recent.faa'
soil_HGT_nonrecent_faa = 'soil_HGT_nonrecent.faa'


# get genome list
Soil_all_taxon_list = set()
Soil_only_taxon_list = set()
for genome_metadata in open(genome_meadata_file):
    genome_metadata_split = genome_metadata.strip().split('\t')
    taxon_id = genome_metadata_split[1]
    niche = genome_metadata_split[2]
    if 'Soil' in niche:
        Soil_all_taxon_list.add(taxon_id)
        if niche == 'Non-Human, Soil':
            Soil_only_taxon_list.add(taxon_id)

print('The number of all soil isolates: %s' % len(Soil_all_taxon_list))


# get func_stats_files for soil isolates
for each_soil_isolate in Soil_all_taxon_list:
    cp_cmd = 'cp func_stats_files_2094_isolates/taxon%s_func_stats.txt func_stats_files_368_soil_isolates/' % each_soil_isolate
    os.system(cp_cmd)


soil_HGT_recent_num = 0
soil_HGT_nonrecent_num = 0
soil_HGT_recent = set()
soil_HGT_nonrecent = set()
soil_HGT_genetic_variation_list = []
n = 0
for each in open(overall_10255_HGTs):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        gene_1_taxon = '_'.join(gene_1.split('_')[:-1])[5:]
        gene_2_taxon = '_'.join(gene_2.split('_')[:-1])[5:]
        genetic_variation = 100 - float(each_split[4])
        direction = each_split[7]

        if (gene_1_taxon in Soil_all_taxon_list) or (gene_2_taxon in Soil_all_taxon_list):

            soil_HGT_genetic_variation_list.append(genetic_variation)

            if (gene_1_taxon in Soil_all_taxon_list) and (gene_2_taxon not in Soil_all_taxon_list):
                if genetic_variation <= 1:
                    soil_HGT_recent.add(gene_1)
                    soil_HGT_recent_num +=1
                else:
                    soil_HGT_nonrecent.add(gene_1)
                    soil_HGT_nonrecent_num += 1
                n += 1

            if (gene_1_taxon not in Soil_all_taxon_list) and (gene_2_taxon in Soil_all_taxon_list):
                if genetic_variation <= 1:
                    soil_HGT_recent.add(gene_2)
                    soil_HGT_recent_num += 1
                else:
                    soil_HGT_nonrecent.add(gene_2)
                    soil_HGT_nonrecent_num += 1
                n += 1

            if (gene_1_taxon in Soil_all_taxon_list) and (gene_2_taxon in Soil_all_taxon_list):
                if genetic_variation <= 1:
                    soil_HGT_recent.add(direction.split('-->')[1])
                    soil_HGT_recent_num += 1
                else:
                    soil_HGT_nonrecent.add(direction.split('-->')[1])
                    soil_HGT_nonrecent_num += 1
                n += 1

print('The number of recent HGT genes from soil isolates: %s' % len(soil_HGT_recent))
print('The number of nonrecent HGT genes from soil isolates: %s' % len(soil_HGT_nonrecent))
print(n)
print(len(soil_HGT_genetic_variation_list))

print('soil_HGT_recent_num: %s' % soil_HGT_recent_num)
print('soil_HGT_nonrecent_num: %s' % soil_HGT_nonrecent_num)

m = 0
for each in soil_HGT_genetic_variation_list:
    if each <= 1:
        m += 1
print('m: %s' % m)

# extract sequences
soil_HGT_recent_faa_handle = open(soil_HGT_recent_faa, 'w')
soil_HGT_nonrecent_faa_handle = open(soil_HGT_nonrecent_faa, 'w')
for gene in SeqIO.parse(overall_10255_HGTs_seq, 'fasta'):

    if gene.id in soil_HGT_recent:
        soil_HGT_recent_faa_handle.write('>%s\n' % gene.id)
        soil_HGT_recent_faa_handle.write('%s\n' % str(gene.seq))

    if gene.id in soil_HGT_nonrecent:
        soil_HGT_nonrecent_faa_handle.write('>%s\n' % gene.id)
        soil_HGT_nonrecent_faa_handle.write('%s\n' % str(gene.seq))

soil_HGT_recent_faa_handle.close()
soil_HGT_nonrecent_faa_handle.close()


# # plot genetic variation
# plt.hist(soil_HGT_genetic_variation_list, bins=30, linewidth=1, alpha=1, color='black', align='mid', rwidth=0.8)
# plt.xticks(range(30))
# plt.show()


#
# # Get plot
# plt.figure(figsize=(10, 5))
# num_bins = range(30)
# plt.hist(soil_HGT_genetic_variation_list, bins=num_bins, alpha=0.6, normed=0, linewidth=0, color='black', rwidth=0.8)
# #plt.title('Genetic variation of identified HGTs')
# #plt.xticks([])
# #plt.xticks(horizontalalignment='center')
# plt.axis('tight')
# bins_labels(range(30), fontsize=12)
#
# plt.xlabel('Genetic divergence (%)')
# plt.ylabel('Number of HGT')
# plt.tight_layout()
# plt.savefig('/Users/songweizhi/Desktop/Figure_8_soil_HGTs_genetic_divergence.png', bbox_inches='tight', dpi=300)
# plt.close()
# plt.clf()
#
#
#
#
#
#
#
#

