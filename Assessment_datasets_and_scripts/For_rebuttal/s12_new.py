import os
from Bio import SeqIO
import shutil
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def bins_labels(bins, **kwargs):
    bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
    #plt.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins, **kwargs)
    plt.xticks(np.arange(min(bins), max(bins), bin_w), bins, **kwargs)
    plt.xlim(bins[0], bins[-1])


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


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


wd = '/Users/songweizhi/Desktop/soil_HGT_377'
os.chdir(wd)


MetaCHIP_hgts = 'Total2094_g664_HGTs_PG_validated.txt'
soil_genome_id = 'taxon_id_Soil368.txt'
MetaCHIP_hgts_nc_seq = 'Total2094_g664_HGTs_PG_nc.fasta'

soil_HGT_seq_folder = 'soil_HGT_seqences_lt97'

soil_HGT_recent_seq = '%s/soil_HGT_recent.faa' % soil_HGT_seq_folder
soil_HGT_nonrecent_seq = '%s/soil_HGT_nonrecent.faa' % soil_HGT_seq_folder
soil_HGT_genetic_variation_plot = 'soil_HGT_genetic_variation.png'

subject_seq = 'Soil_HGT_sequences_by_id.fasta'
paired_16s_blast_out = 'output_min700bp.tab'


blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'



# get 16s identity dict
iden_dict = {}
for each in open(paired_16s_blast_out):
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


soil_HGT_iden_list = []
soil_HGT_recent = set()
soil_HGT_nonrecent = set()
soil_HGT_recent_num = 0
soil_HGT_nonrecent_num = 0
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

        # only work on HGTs with 16s less than 97
        if iden_16s < 97:

            if (recipient_genome in soil_genome_id_list) and (donor_genome in soil_genome_id_list):
                gene_to_keep = ''
                if gene_1_genome == recipient_genome:
                    gene_to_keep = gene_1
                if gene_2_genome == recipient_genome:
                    gene_to_keep = gene_2

                if identity >= 99:
                    soil_HGT_recent.add(gene_to_keep)
                    soil_HGT_recent_num += 1
                else:
                    soil_HGT_nonrecent.add(gene_to_keep)
                    soil_HGT_nonrecent_num += 1

                soil_HGT_iden_list.append(identity)

            elif (recipient_genome in soil_genome_id_list) and (donor_genome not in soil_genome_id_list):
                gene_to_keep = ''
                if gene_1_genome == recipient_genome:
                    gene_to_keep = gene_1
                if gene_2_genome == recipient_genome:
                    gene_to_keep = gene_2

                if identity >= 99:
                    soil_HGT_recent.add(gene_to_keep)
                    soil_HGT_recent_num += 1
                else:
                    soil_HGT_nonrecent.add(gene_to_keep)
                    soil_HGT_nonrecent_num += 1

                soil_HGT_iden_list.append(identity)

            elif (recipient_genome not in soil_genome_id_list) and (donor_genome in soil_genome_id_list):
                gene_to_keep = ''
                if gene_1_genome == donor_genome:
                    gene_to_keep = gene_1
                if gene_2_genome == donor_genome:
                    gene_to_keep = gene_2

                if identity >= 99:
                    soil_HGT_recent.add(gene_to_keep)
                    soil_HGT_recent_num += 1
                else:
                    soil_HGT_nonrecent.add(gene_to_keep)
                    soil_HGT_nonrecent_num += 1

                soil_HGT_iden_list.append(identity)


print('The number of recent HGTs involving soil genome: %s' % soil_HGT_recent_num)
print('The number of nonrecent HGTs involving soil genome: %s' % soil_HGT_nonrecent_num)


force_create_folder(soil_HGT_seq_folder)

# extract HGT sequences
soil_HGT_recent_seq_handle = open(soil_HGT_recent_seq, 'w')
soil_HGT_nonrecent_seq_handle = open(soil_HGT_nonrecent_seq, 'w')
hgts_seq_involving_soil_handle = open('hgts_seq_involving_soil.fasta', 'w')
for each_seq in SeqIO.parse(MetaCHIP_hgts_nc_seq, 'fasta'):

    if str(each_seq.id) in soil_HGT_recent:
        soil_HGT_recent_seq_handle.write('>%s\n' % str(each_seq.id))
        soil_HGT_recent_seq_handle.write('%s\n' % str(each_seq.seq))

    if str(each_seq.id) in soil_HGT_nonrecent:
        soil_HGT_nonrecent_seq_handle.write('>%s\n' % str(each_seq.id))
        soil_HGT_nonrecent_seq_handle.write('%s\n' % str(each_seq.seq))

soil_HGT_recent_seq_handle.close()
soil_HGT_nonrecent_seq_handle.close()


# COG annotation of soil HGTs
COG_cmd_recent_hgt =    'MyBioTools COG_annot -i %s -x faa -m N -p %s -t 2 -DB_dir /Users/songweizhi/COG_DB' % (soil_HGT_seq_folder, soil_HGT_seq_folder)
#os.system(COG_cmd_recent_hgt)


######################################## plot identity list ########################################

# transfer identity to genetic variation
soil_HGT_genetic_variation_list = []
for each in soil_HGT_iden_list:
    soil_HGT_genetic_variation_list.append(100 -each)


# Get plot
num_bins = range(30)
plt.hist(soil_HGT_genetic_variation_list, bins=num_bins, alpha=0.8, normed=0, linewidth=0, color='black', rwidth=0.8)
#plt.title('Genetic variation of identified HGTs')
#plt.xticks([])
#plt.xticks(horizontalalignment='center')
plt.axis('tight')
bins_labels(range(30), fontsize=15)

plt.xlabel('Genetic divergence (%)')
plt.ylabel('Number of HGT')
plt.tight_layout()
plt.savefig(soil_HGT_genetic_variation_plot, bbox_inches='tight', dpi=300)
plt.close()
plt.clf()













# # run blast between genes transferred and pubished sequence fragments
# blast_output_hgts_involving_soil = 'blast_output_hgts_involving_soil.tab'
# blast_command_involving_soil = 'blastn -query %s -subject %s -out %s %s' % ('hgts_seq_involving_soil.fasta', subject_seq, blast_output_hgts_involving_soil, blast_parameters)
# os.system(blast_command_involving_soil)
#
# # parse blast results for recent HGTs
# overlapped_predictions_involving_soil = get_overlapped_predictions(blast_output_hgts_involving_soil)
# #print(len(overlapped_predictions_involving_soil))
# print('The number of overlapped predictions: %s' % len(overlapped_predictions_involving_soil))
#
#
# # get sequence length dict
# gene_length_dict = {}
# for each_seq in SeqIO.parse(MetaCHIP_hgts_nc_seq, 'fasta'):
#     gene_length_dict[str(each_seq.id)] = len(str(each_seq.seq))
#

# for_flk_check = set()
# n = 0
# shorter500 = 0
# for each_gene in hgts_involving_soil:
#     if each_gene not in overlapped_predictions_involving_soil:
#         n += 1
#         #print('%s\t%s' % (each_gene, gene_length_dict[each_gene]))
#         if gene_length_dict[each_gene] <= 500:
#             shorter500 += 1
#             for_flk_check.add(each_gene)
#
#
