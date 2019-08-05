import os
from Bio import SeqIO


# input faa file
genome_metadata = '/Users/songweizhi/Desktop/nature_rebuttal/nature10571_genome_metadata_2094.txt'
HGT_BM_faa_file = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_BM_aa.fasta'
HGT_PG_faa_file = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_aa.fasta'


human_associated_taxons = set()
for genome in open(genome_metadata):
    if not genome.startswith('Genome_Name'):
        genome_split = genome.strip().split('\t')
        taxon = genome_split[1]
        enviro = genome_split[2]
        if ('Non-Human' not in enviro) and ('Soil' not in enviro) and ('Marine' not in enviro) and ('Airways' not in enviro):
            human_associated_taxons.add(taxon)
print(len(human_associated_taxons))


for each_taxon in human_associated_taxons:
    cp_cmd = 'cp /Users/songweizhi/Desktop/func_stats_files/taxon%s_func_stats.txt /Users/songweizhi/Desktop/human_associated_func_stats_files/' % each_taxon
    #os.system(cp_cmd)


# output human associated faa file
Human_associated_HGT_BM_faa_file = '/Users/songweizhi/Desktop/Human_associated_Total2094_g664_HGTs_BM_aa.fasta'
Human_associated_HGT_PG_faa_file = '/Users/songweizhi/Desktop/Human_associated_Total2094_g664_HGTs_PG_aa.fasta'


Human_associated_HGT_BM_faa_file_handle = open(Human_associated_HGT_BM_faa_file, 'w')
for BM_HGT in SeqIO.parse(HGT_BM_faa_file, 'fasta'):
    BM_HGT_id = BM_HGT.id
    BM_HGT_taxon = BM_HGT_id.split('_')[0][5:]
    if BM_HGT_taxon in human_associated_taxons:
        Human_associated_HGT_BM_faa_file_handle.write('>%s\n' % BM_HGT_id)
        Human_associated_HGT_BM_faa_file_handle.write('%s\n' % str(BM_HGT.seq))
Human_associated_HGT_BM_faa_file_handle.close()


Human_associated_HGT_PG_faa_file_handle = open(Human_associated_HGT_PG_faa_file, 'w')
for PG_HGT in SeqIO.parse(HGT_PG_faa_file, 'fasta'):
    PG_HGT_id = PG_HGT.id
    PG_HGT_taxon = PG_HGT_id.split('_')[0][5:]
    if PG_HGT_taxon in human_associated_taxons:
        Human_associated_HGT_PG_faa_file_handle.write('>%s\n' % PG_HGT_id)
        Human_associated_HGT_PG_faa_file_handle.write('%s\n' % str(PG_HGT.seq))
Human_associated_HGT_PG_faa_file_handle.close()
