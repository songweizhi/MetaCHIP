import os
from Bio import SeqIO
import matplotlib.pyplot as plt


wd = '/Users/songweizhi/Desktop/cor'

os.chdir(wd)


bin_list = []
bin_hgt_num_dict = {}
for each in open('GoodBins_0.5_0.05_grouping_o34.txt'):
    bin_id = each.strip().split(',')[1]
    bin_list.append(bin_id)
    bin_hgt_num_dict[bin_id] = 0


bin_contamination_dict = {}
for each in open('bin_conta.txt'):
    each_split = each.strip().split('\t')
    bin_contamination_dict[each_split[0]] = float(each_split[1])

n = 0
HGT_results = 'HGT_candidates_PG_validated.txt'
HGT_results = 'HGT_candidates_PG.txt'

#HGT_results = 'HGT_candidates_BM.txt'
for each in open(HGT_results):

    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        iden = float(each_split[4])
        pg_validated = each_split[6]
        gene_1 = each_split[0]
        genome_1 = '_'.join(gene_1.split('_')[:-1])
        gene_2 = each_split[1]
        genome_2 = '_'.join(gene_2.split('_')[:-1])

        if (iden >= 100) and (pg_validated != 'N/A'):
            print(each.strip())
            gene_1_ctg_len = 0
            gene_2_ctg_len = 0
            genome_1_gbk = 'GoodBins_0.5_0.05_gbk_files/%s.gbk' % genome_1
            for genome_1_record in SeqIO.parse(genome_1_gbk, 'genbank'):
                for gene in genome_1_record.features:
                    if 'locus_tag' in gene.qualifiers:
                        if gene_1 in gene.qualifiers["locus_tag"]:
                            gene_1_ctg_len = len(genome_1_record.seq)

            genome_2_gbk = 'GoodBins_0.5_0.05_gbk_files/%s.gbk' % genome_2
            for genome_2_record in SeqIO.parse(genome_2_gbk, 'genbank'):
                for gene in genome_2_record.features:
                    if 'locus_tag' in gene.qualifiers:
                        if gene_2 in gene.qualifiers["locus_tag"]:
                            gene_2_ctg_len = len(genome_2_record.seq)

            if gene_1_ctg_len > gene_2_ctg_len:
                bin_hgt_num_dict[genome_2] += 1

            if gene_2_ctg_len > gene_1_ctg_len:
                bin_hgt_num_dict[genome_1] += 1

            print(n)
            n += 1


contamination_list = []
HGT_num_list = []
for each in bin_list:
    contamination_list.append(bin_contamination_dict[each])
    HGT_num_list.append(bin_hgt_num_dict[each])


plt.scatter(contamination_list, HGT_num_list, alpha=0.05, linewidth=0)
plt.xlabel("Contamination (%)")
plt.ylabel("HGT number")
plt.show()

