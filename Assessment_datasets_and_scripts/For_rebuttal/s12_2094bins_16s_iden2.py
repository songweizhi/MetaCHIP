import os
from Bio import SeqIO


def get_transposase_related_COGs(cog_stats_file):
    transposase_related_COGs = set()
    for each_cog in open(cog_stats_file):
        each_cog_split = each_cog.strip().split('\t')
        cog_id = each_cog_split[0]
        cog_description = each_cog_split[1]
        if 'Transposase' in cog_description:
            transposase_related_COGs.add(cog_id)

    return transposase_related_COGs


def get_transposase_related_genes(protein_id_cog_file, Transposase_related_COGs):

    Transposase_related_genes = set()
    for each_protein in open(protein_id_cog_file):
        each_protein_split = each_protein.strip().split('\t')
        protein_id = each_protein_split[0]
        cog_id = each_protein_split[1]
        if cog_id in Transposase_related_COGs:
            Transposase_related_genes.add(protein_id)

    return Transposase_related_genes


def get_non_transposase_hgts(hgts_genes, Transposase_related_genes):

    hgts_no_transposase = set()
    for each_hgt_gene in hgts_genes:
        if each_hgt_gene not in Transposase_related_genes:
            hgts_no_transposase.add(each_hgt_gene)

    return hgts_no_transposase


def get_gene_sequences(seq_file_in, gene_id_list, seq_file_out):
    seq_file_out_handle = open(seq_file_out, 'w')
    for each_seq in SeqIO.parse(seq_file_in, 'fasta'):
        if str(each_seq.id) in gene_id_list:
            seq_file_out_handle.write('>%s\n' % str(each_seq.id))
            seq_file_out_handle.write('%s\n' % str(each_seq.seq))
    seq_file_out_handle.close()


def get_overlapped_predictions(blast_output):
    overlapped_predictions = set()
    for blast_hit in open(blast_output):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        identity = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_len = int(blast_hit_split[12])
        query_cov = align_len / query_len
        if (identity == 100) and (query_cov == 1):
            overlapped_predictions.add(query)

    return overlapped_predictions


wd = '/Users/songweizhi/Desktop/test'
os.chdir(wd)


MetaCHIP_hgts = 'Total2094_g664_HGTs_PG_validated.txt'
soil_genome_id = '/Users/songweizhi/Desktop/nature_rebuttal/taxon_id_Soil368.txt'
MetaCHIP_hgts_nc_seq = 'Total2094_g664_HGTs_PG_nc.fasta'
subject_seq = '/Users/songweizhi/Desktop/nature_rebuttal/Soil_HGT_sequences_by_id.fasta'
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'



# get 16s identity dict
iden_dict = {}
for each in open('/Users/songweizhi/Desktop/output_min700bp.tab'):
    each_split = each.strip().split('\t')
    genome_1 = each_split[0].split('_')[0]
    genome_2 = each_split[1].split('_')[0]
    key = '%s_%s' % (genome_1, genome_2)
    aln = int(each_split[3])
    iden = float(each_split[2])
    if aln >= 700:
        if key not in iden_dict:
            iden_dict[key] = iden
        else:
            if iden > iden_dict[key]:
                iden_dict[key] = iden


soil_genome_id_list = set()
for soil_genome in open(soil_genome_id):
    soil_genome_id_list.add(soil_genome.strip())


hgts_involving_soil = set()

recent_num = 0
total_num = 0
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
        key_value = '%s_%s' % (gene_1_genome, gene_2_genome)
        iden_16s = 0
        if key_value in iden_dict:
            iden_16s = iden_dict[key_value]

        if (identity > 99):
            recent_num += 1

        total_num += 1

print(total_num)
print(recent_num)

