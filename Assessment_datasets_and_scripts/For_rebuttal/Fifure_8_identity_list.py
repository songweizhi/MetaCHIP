import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def bins_labels(bins, **kwargs):
    bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
    plt.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins, **kwargs)
    plt.xlim(bins[0], bins[-1])


# get identity list
identity_list = []
for hgt_pcofg in open('/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_validated.txt'):
    if not hgt_pcofg.startswith('Gene_1'):
        hgt_pcofg_split = hgt_pcofg.strip().split('\t')
        concatenated_genes = '%s___%s' % (hgt_pcofg_split[0], hgt_pcofg_split[1])
        identity = float(hgt_pcofg_split[4])
        identity_list.append(identity)


# transfer identity to genetic variation
genetic_variation_list = []
for each in identity_list:
    genetic_variation_list.append(100 -each)


# Get plot
num_bins = range(30)
plt.hist(genetic_variation_list, bins=num_bins, alpha=0.6, normed=0, linewidth=0, color='black', rwidth=0.8)
#plt.title('Genetic variation of identified HGTs')
#plt.xticks([])
#plt.xticks(horizontalalignment='center')
plt.axis('tight')
bins_labels(range(30), fontsize=12)

plt.xlabel('Genetic variation (%)')
plt.ylabel('Number of HGT')
plt.tight_layout()
plt.savefig('/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_validated.png', bbox_inches='tight', dpi=300)
plt.close()
plt.clf()



