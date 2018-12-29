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
MetaCHIP_hgts_aa_seq = 'Total2094_g664_HGTs_PG_aa.fasta'
MetaCHIP_hgts_nc_seq = 'Total2094_g664_HGTs_PG_nc.fasta'
subject_seq = '/Users/songweizhi/Desktop/nature_rebuttal/nature10571_hgt_seqs/hgt.fst'
#subject_seq = '/Users/songweizhi/Desktop/nature_rebuttal/Soil_HGT_sequences_by_id.fasta'
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'


# get 16s identity dict
iden_dict = {}
for each in open('/Users/songweizhi/Desktop/output.tab'):
    each_split = each.strip().split('\t')
    genome_1 = each_split[0].split('_')[0]
    genome_2 = each_split[1].split('_')[0]
    key = '%s_%s' % (genome_1, genome_2)
    aln = int(each_split[3])
    iden = float(each_split[2])
    if aln >= 1300:
        if key not in iden_dict:
            iden_dict[key] = iden
        else:
            if iden > iden_dict[key]:
                iden_dict[key] = iden


soil_genome_id_list = set()
for soil_genome in open(soil_genome_id):
    soil_genome_id_list.add(soil_genome.strip())


hgts_involving_soil = set()
# hgts_to_soil = set()
# hgts_from_soil = set()
# hgts_within_soil = set()
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


        if (identity >= 99) and (iden_16s < 97):

            # get genes transferred within soil genomes
            if (recipient_genome in soil_genome_id_list) and (donor_genome in soil_genome_id_list):
                if gene_1_genome == recipient_genome:
                    #hgts_within_soil.add(gene_1)
                    hgts_involving_soil.add(gene_1)
                if gene_2_genome == recipient_genome:
                    #hgts_within_soil.add(gene_2)
                    hgts_involving_soil.add(gene_2)
            else:
                # get genes transferred to soil genomes
                if (recipient_genome in soil_genome_id_list) and (donor_genome not in soil_genome_id_list):
                    if gene_1_genome == recipient_genome:
                        #hgts_to_soil.add(gene_1)
                        hgts_involving_soil.add(gene_1)
                    if gene_2_genome == recipient_genome:
                        #hgts_to_soil.add(gene_2)
                        hgts_involving_soil.add(gene_2)

                # get genes transferred from soil genomes
                if (donor_genome in soil_genome_id_list) and (recipient_genome not in soil_genome_id_list):
                    if gene_1_genome == donor_genome:
                        #hgts_from_soil.add(gene_1)
                        hgts_involving_soil.add(gene_1)
                    if gene_2_genome == donor_genome:
                        #hgts_from_soil.add(gene_2)
                        hgts_involving_soil.add(gene_2)

#print('The number of soil genomes: %s' % len(soil_genome_id_list))
print('The number of recent HGTs involving soil genome: %s' % len(hgts_involving_soil))
# print('The number of recent HGTs transferred to soil genome: %s' % len(hgts_to_soil))
# print('The number of recent HGTs transferred from soil genome: %s' % len(hgts_from_soil))
# print('The number of recent HGTs transferred within soil genome: %s' % len(hgts_within_soil))
#

# extract HGT sequences
hgts_seq_involving_soil_handle = open('hgts_seq_involving_soil.fasta', 'w')
# hgts_seq_to_soil_handle = open('hgts_seq_to_soil.fasta', 'w')
# hgts_seq_from_soil_handle = open('hgts_seq_from_soil.fasta', 'w')
# hgts_seq_within_soil_handle = open('hgts_seq_within_soil.fasta', 'w')
for each_seq in SeqIO.parse(MetaCHIP_hgts_nc_seq, 'fasta'):
    if str(each_seq.id) in hgts_involving_soil:
        hgts_seq_involving_soil_handle.write('>%s\n' % str(each_seq.id))
        hgts_seq_involving_soil_handle.write('%s\n' % str(each_seq.seq))
    # if str(each_seq.id) in hgts_to_soil:
    #     hgts_seq_to_soil_handle.write('>%s\n' % str(each_seq.id))
    #     hgts_seq_to_soil_handle.write('%s\n' % str(each_seq.seq))
    # if str(each_seq.id) in hgts_from_soil:
    #     hgts_seq_from_soil_handle.write('>%s\n' % str(each_seq.id))
    #     hgts_seq_from_soil_handle.write('%s\n' % str(each_seq.seq))
    # if str(each_seq.id) in hgts_within_soil:
    #     hgts_seq_within_soil_handle.write('>%s\n' % str(each_seq.id))
    #     hgts_seq_within_soil_handle.write('%s\n' % str(each_seq.seq))
hgts_seq_involving_soil_handle.close()
# hgts_seq_to_soil_handle.close()
# hgts_seq_from_soil_handle.close()
# hgts_seq_within_soil_handle.close()


# get Transposase related COG id
hgts_involving_soil_cog_stats_file = 'hgts_seq_involving_soil_cog_stats.txt'
# hgts_to_soil_cog_stats_file =        'hgts_seq_to_soil_cog_stats.txt'
# hgts_from_soil_cog_stats_file =      'hgts_seq_from_soil_cog_stats.txt'
# hgts_within_soil_cog_stats_file =    'hgts_seq_within_soil_cog_stats.txt'
Transposase_related_COGs_involving_soil = get_transposase_related_COGs(hgts_involving_soil_cog_stats_file)
# Transposase_related_COGs_to_soil =        get_transposase_related_COGs(hgts_to_soil_cog_stats_file)
# Transposase_related_COGs_from_soil =      get_transposase_related_COGs(hgts_from_soil_cog_stats_file)
# Transposase_related_COGs_within_soil =    get_transposase_related_COGs(hgts_within_soil_cog_stats_file)


# get Transposase related gene id
protein_id_cog_file_involving_soil = 'protein-id_cog_involving_soil.txt'
# protein_id_cog_file_to_soil =        'protein-id_cog_to_soil.txt'
# protein_id_cog_file_from_soil =      'protein-id_cog_from_soil.txt'
# protein_id_cog_file_within_soil =    'protein-id_cog_within_soil.txt'
Transposase_related_genes_involving_soil = get_transposase_related_genes(protein_id_cog_file_involving_soil, Transposase_related_COGs_involving_soil)
# Transposase_related_genes_to_soil =        get_transposase_related_genes(protein_id_cog_file_to_soil, Transposase_related_COGs_to_soil)
# Transposase_related_genes_from_soil =      get_transposase_related_genes(protein_id_cog_file_from_soil, Transposase_related_COGs_from_soil)
# Transposase_related_genes_within_soil =    get_transposase_related_genes(protein_id_cog_file_within_soil, Transposase_related_COGs_within_soil)


# get genes transferred (no transposase)
hgts_involving_soil_no_transposase = get_non_transposase_hgts(hgts_involving_soil, Transposase_related_genes_involving_soil)
# hgts_to_soil_no_transposase =        get_non_transposase_hgts(hgts_to_soil, Transposase_related_genes_to_soil)
# hgts_from_soil_no_transposase =      get_non_transposase_hgts(hgts_from_soil, Transposase_related_genes_from_soil)
# hgts_within_soil_no_transposase =    get_non_transposase_hgts(hgts_within_soil, Transposase_related_genes_within_soil)
# print('The number of recent HGTs involving soil genome (no_transposase): %s' % len(hgts_involving_soil_no_transposase))


# get sequence of transferred genes (no transposase)
hgts_involving_soil_no_transposase_seq_file = 'hgts_seq_involving_soil_no_transposase.fasta'
# hgts_to_soil_no_transposase_seq_file = 'hgts_seq_to_soil_no_transposase.fasta'
# hgts_from_soil_no_transposase_seq_file = 'hgts_seq_from_soil_no_transposase.fasta'
# hgts_within_soil_no_transposase_seq_file = 'hgts_seq_within_soil_no_transposase.fasta'
get_gene_sequences(MetaCHIP_hgts_nc_seq, hgts_involving_soil_no_transposase, hgts_involving_soil_no_transposase_seq_file)
# get_gene_sequences(MetaCHIP_hgts_nc_seq, hgts_to_soil_no_transposase, hgts_to_soil_no_transposase_seq_file)
# get_gene_sequences(MetaCHIP_hgts_nc_seq, hgts_from_soil_no_transposase, hgts_from_soil_no_transposase_seq_file)
# get_gene_sequences(MetaCHIP_hgts_nc_seq, hgts_within_soil_no_transposase, hgts_within_soil_no_transposase_seq_file)
#

# run blast between genes transferred (no transposase) and pubished soil sequence fragments
blast_output_hgts_involving_soil =  'blast_output_hgts_involving_soil.tab'
# blast_output_hgts_to_soil =         'blast_output_hgts_to_soil.tab'
# blast_output_hgts_from_soil =       'blast_output_hgts_from_soil.tab'
# blast_output_hgts_within_soil =     'blast_output_hgts_within_soil.tab'

blast_command_involving_soil = 'blastn -query %s -subject %s -out %s %s' % ('hgts_seq_involving_soil.fasta', subject_seq, blast_output_hgts_involving_soil, blast_parameters)
# blast_command_to_soil =        'blastn -query %s -subject %s -out %s %s' % ('hgts_seq_to_soil.fasta',        subject_seq, blast_output_hgts_to_soil,        blast_parameters)
# blast_command_from_soil =      'blastn -query %s -subject %s -out %s %s' % ('hgts_seq_from_soil.fasta',      subject_seq, blast_output_hgts_from_soil,      blast_parameters)
# blast_command_withim_soil =    'blastn -query %s -subject %s -out %s %s' % ('hgts_seq_within_soil.fasta',    subject_seq, blast_output_hgts_within_soil,    blast_parameters)

# blast_command_involving_soil = 'blastn -query %s -subject %s -out %s %s' % (hgts_involving_soil_no_transposase_seq_file, subject_seq, blast_output_hgts_involving_soil, blast_parameters)
# blast_command_to_soil =        'blastn -query %s -subject %s -out %s %s' % (hgts_to_soil_no_transposase_seq_file,        subject_seq, blast_output_hgts_to_soil,        blast_parameters)
# blast_command_from_soil =      'blastn -query %s -subject %s -out %s %s' % (hgts_from_soil_no_transposase_seq_file,      subject_seq, blast_output_hgts_from_soil,      blast_parameters)
# blast_command_withim_soil =    'blastn -query %s -subject %s -out %s %s' % (hgts_within_soil_no_transposase_seq_file,    subject_seq, blast_output_hgts_within_soil,    blast_parameters)

#os.system(blast_command_involving_soil)
#os.system(blast_command_to_soil)
#os.system(blast_command_from_soil)
#os.system(blast_command_withim_soil)


# parse blast results for recent HGTs
overlapped_predictions_involving_soil = get_overlapped_predictions(blast_output_hgts_involving_soil)
# overlapped_predictions_to_soil        = get_overlapped_predictions(blast_output_hgts_to_soil)
# overlapped_predictions_from_soil      = get_overlapped_predictions(blast_output_hgts_from_soil)
# overlapped_predictions_within_soil    = get_overlapped_predictions(blast_output_hgts_within_soil)

print(len(overlapped_predictions_involving_soil))
# print(len(overlapped_predictions_to_soil))
# print(len(overlapped_predictions_from_soil))
# print(len(overlapped_predictions_within_soil))


# get sequence length dict
gene_length_dict = {}
for each_seq in SeqIO.parse(MetaCHIP_hgts_nc_seq, 'fasta'):
    gene_length_dict[str(each_seq.id)] = len(str(each_seq.seq))


m = 0
n = 0
for each_gene in hgts_involving_soil:

    if each_gene in overlapped_predictions_involving_soil:

        #print(each_gene)
        m += 1
    if each_gene not in overlapped_predictions_involving_soil:

        #print(each_gene)
        n += 1

print(m)
print(n)
#
# n = 0
# for each_hgt in open(MetaCHIP_hgts):
#     if not each_hgt.startswith('Gene_1\tGene_2'):
#         each_hgt_split = each_hgt.strip().split('\t')
#         gene_1 = each_hgt_split[0]
#         gene_2 = each_hgt_split[1]
#         identity = float(each_hgt_split[4])
#         direction = each_hgt_split[7]
#         gene_1_genome = '_'.join(gene_1.split('_')[:-1])
#         gene_2_genome = '_'.join(gene_2.split('_')[:-1])
#         donor_genome = direction.split('-->')[0]
#         recipient_genome = direction.split('-->')[1]
#         if identity >= 99:
#
#             if (gene_1_genome in soil_genome_id_list) or (gene_2_genome in soil_genome_id_list):
#
#                 if (gene_1 in overlapped_predictions_involving_soil) or (gene_2 in overlapped_predictions_involving_soil):
#                     n += 1
#
#                 #print(each_hgt)
#
#
#                 # key_value = '%s_%s' % (gene_1_genome, gene_2_genome)
#                 # iden_16s = 0
#                 # if key_value in iden_dict:
#                 #     iden_16s = iden_dict[key_value]
#                 #
#                 # if (gene_1 in overlapped_predictions_involving_soil) or (gene_2 in overlapped_predictions_involving_soil):
#                 #     pass
#                 #     print('%s\t%s\t%s' % (each_hgt.strip(), 'yes', iden_16s))
#                 # else:
#                 #
#                 #
#                 #     if iden_16s >=97:
#                 #         #print('%s\t%s\t%s' % (each_hgt.strip(), 'no', iden_16s))
#                 #         n += 1
#
# print(n)
#
#
