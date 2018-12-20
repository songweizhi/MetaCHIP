import os
from Bio import SeqIO


MetaCHIP_hgts = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_validated.txt'
#MetaCHIP_hgts = '/Users/songweizhi/Desktop/Total2094_f285_HGTs_PG_validated.txt'
MetaCHIP_hgts_aa_seq = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_aa.fasta'
MetaCHIP_hgts_nc_seq = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_nc.fasta'

soil_genome_id = '/Users/songweizhi/Desktop/nature_rebuttal/taxon_id_Soil368.txt'

soil_genome_id_list = set()
for soil_genome in open(soil_genome_id):
    soil_genome_id_list.add(soil_genome.strip())


soil_hgts = set()
for each_hgt in open(MetaCHIP_hgts):

    if not each_hgt.startswith('Gene_1\tGene_2'):

        each_hgt_split = each_hgt.strip().split('\t')
        gene_1 = each_hgt_split[0]
        gene_2 = each_hgt_split[1]
        identity = float(each_hgt_split[4])
        direction = each_hgt_split[7]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        donor_genome = direction.split('-->')[0]
        recipient_genome = direction.split('-->')[1]

        if identity >= 99:

            if gene_1_genome in soil_genome_id_list:
                soil_hgts.add(gene_1)

            if gene_2_genome in soil_genome_id_list:
                soil_hgts.add(gene_2)


print('The number of soil genomes: %s' % len(soil_genome_id_list))
print('The number of recent HGTs from soil genome: %s' % len(soil_hgts))


MetaCHIP_hgts_aa_seq_soil = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_aa_soil623.fasta'
MetaCHIP_hgts_aa_seq_soil_handle = open(MetaCHIP_hgts_aa_seq_soil, 'w')
for each_aa in SeqIO.parse(MetaCHIP_hgts_aa_seq, 'fasta'):
    if str(each_aa.id) in soil_hgts:
        MetaCHIP_hgts_aa_seq_soil_handle.write('>%s\n' % str(each_aa.id))
        MetaCHIP_hgts_aa_seq_soil_handle.write('%s\n' % str(each_aa.seq))
MetaCHIP_hgts_aa_seq_soil_handle.close()


# get Transposase related COG id
cog_stats_file = '/Users/songweizhi/Desktop/soil623_cog_stats.txt'
Transposase_related_COGs = set()
for each_cog in open(cog_stats_file):
    each_cog_split = each_cog.strip().split('\t')
    cog_id = each_cog_split[0]
    cog_description = each_cog_split[1]
    if 'Transposase' in cog_description:
        Transposase_related_COGs.add(cog_id)


# get Transposase related HGTs
Transposase_related_genes = set()
protein_id_cog_file = '/Users/songweizhi/Desktop/soil623_protein-id_cog.txt'
for each_protein in open(protein_id_cog_file):
    each_protein_split = each_protein.strip().split('\t')
    protein_id = each_protein_split[0]
    cog_id = each_protein_split[1]
    if cog_id in Transposase_related_COGs:
        Transposase_related_genes.add(protein_id)


# get genes transferred to soil genomes (no transposase)
hgts_to_soil_genome_no_transposase = set()
for each_hgt_gene in soil_hgts:
    if each_hgt_gene not in Transposase_related_genes:
        hgts_to_soil_genome_no_transposase.add(each_hgt_gene)
print('The number of recent HGTs transferred to soil genome (transposase): %s' % len(Transposase_related_genes))
print('The number of recent HGTs transferred to soil genome (no transposase): %s' % len(hgts_to_soil_genome_no_transposase))


# get sequence of genes transferred to soil genomes (no transposase)
MetaCHIP_hgts_aa_seq_soil_no_transposase = '/Users/songweizhi/Desktop/Total2094_g664_HGTs_PG_aa_to_soil_no_transposase.fasta'
MetaCHIP_hgts_aa_seq_soil_no_transposase_handle = open(MetaCHIP_hgts_aa_seq_soil_no_transposase, 'w')
for each_aa in SeqIO.parse(MetaCHIP_hgts_nc_seq, 'fasta'):
    if str(each_aa.id) in hgts_to_soil_genome_no_transposase:
        MetaCHIP_hgts_aa_seq_soil_no_transposase_handle.write('>%s\n' % str(each_aa.id))
        MetaCHIP_hgts_aa_seq_soil_no_transposase_handle.write('%s\n' % str(each_aa.seq))
MetaCHIP_hgts_aa_seq_soil_no_transposase_handle.close()


# run blast between genes transferred to soil genomes (no transposase) and pubished soil sequence fragments
recent_soil_hgt_blast_output = '/Users/songweizhi/Desktop/HGTs_to_soil_no_transposase.tab'
#subject_seq = '/Users/songweizhi/Desktop/nature_rebuttal/Soil_HGT_sequences_by_id.fasta'
subject_seq = '/Users/songweizhi/Desktop/nature_rebuttal/nature10571_hgt_seqs/hgt.fst'
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
recent_hgt_blast_command = 'blastn -query %s -subject %s -out %s %s' % (MetaCHIP_hgts_aa_seq_soil_no_transposase, subject_seq, recent_soil_hgt_blast_output, blast_parameters)
os.system(recent_hgt_blast_command)


# parse blast results for recent HGTs (should blast genes from recipients only, current both donor and recipient genes)
validated_recent_hgt = set()
for blast_hit in open(recent_soil_hgt_blast_output):
    blast_hit_split = blast_hit.strip().split('\t')
    query = blast_hit_split[0]
    identity = float(blast_hit_split[2])
    align_len = int(blast_hit_split[3])
    query_len = int(blast_hit_split[12])
    subject_len = int(blast_hit_split[13])
    query_cov = align_len / query_len
    subject_cov = align_len / subject_len
    if (identity == 100) and (query_cov == 1):
        validated_recent_hgt.add(query)


print(len(validated_recent_hgt))



#
#
#
#




















