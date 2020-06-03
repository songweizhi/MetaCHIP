import os
from Bio import SeqIO


def check_list_elements_exist(elements_to_check, query_object):

    elements_exist = 0
    for each_element in elements_to_check:
        if each_element in query_object:
            elements_exist = 1

    return elements_exist


# file in
HGT_seq_file =                      '/Users/songweizhi/Desktop/rrr/SpongeMAGs_dRep99_pcofg_detected_HGTs_recipient_genes.faa'
prodigal_output_folder =            '/Users/songweizhi/Desktop/rrr/prodigal_output'
recipient_gene_same_ctg_gene_seqs = '/Users/songweizhi/Desktop/rrr/recipient_gene_same_ctg_genes.fasta'
KEGG_annotation_results =           '/Users/songweizhi/Desktop/rrr/recipient_gene_same_ctg_genes_KEGG_wd/recipient_gene_same_ctg_genes_ko_assignment_ABCD.txt'

key_word_list_C = ['DNA repair and recombination proteins', 'Homologous recombination']
key_word_list_D = ['Phage',         'phage',
                   'Transposase',   'transposase', 'transposon',
                   'Integrase',     'integrase',
                   'Recombinase',   'recombinase',
                   'Plasmid',       'plasmid',
                   'Recombination', 'recombination',
                   'insertion element']

'''
KEGG: 
exonuclease
'''

# db file
COG_db_dir =                        '/Users/songweizhi/DB/COG2014'
KEGG_db_dir =                       '/Users/songweizhi/DB/KEGG_2016-09-26'




# get recipient gene list
recipient_gene_list = []
for gene in SeqIO.parse(HGT_seq_file, 'fasta'):
    recipient_gene_list.append(gene.id)

# get contigs which the identified HGTs sit in
recipient_gene_same_ctg_gene_seqs_handle = open(recipient_gene_same_ctg_gene_seqs, 'w')
recipient_gene_ctg_dict = {}
recipient_gene_same_ctg_genes_dict = {}
for recipient_gene in recipient_gene_list:
    recipient_genome = '_'.join(recipient_gene.split('_')[:-1])
    recipient_genome_gbk = '%s/%s.gbk' % (prodigal_output_folder, recipient_genome)
    for record in SeqIO.parse(recipient_genome_gbk, 'genbank'):
        ctg_id = record.id
        ctg_coded_gene_list = []
        for gene in record.features:
            if (gene.type != 'source') and ('locus_tag' in gene.qualifiers):
                gene_id = gene.qualifiers['locus_tag'][0]
                ctg_coded_gene_list.append(gene_id)

        if recipient_gene in ctg_coded_gene_list:
            for gene in record.features:
                if (gene.type != 'source') and ('locus_tag' in gene.qualifiers):

                    gene_id = gene.qualifiers['locus_tag'][0]
                    gene_seq = gene.qualifiers['translation'][0]
                    recipient_gene_same_ctg_gene_seqs_handle.write('>%s\n' % gene_id)
                    recipient_gene_same_ctg_gene_seqs_handle.write('%s\n' % gene_seq)

            recipient_gene_ctg_dict[recipient_gene] = ctg_id
            recipient_gene_same_ctg_genes_dict[recipient_gene] = ctg_coded_gene_list

recipient_gene_same_ctg_gene_seqs_handle.close()

# annotate flk genes
BioSAK_COG_cmd = 'BioSAK COG2014 -m P -t 4 -db_dir %s -i %s -diamond' % (COG_db_dir, recipient_gene_same_ctg_gene_seqs)
BioSAK_KEGG_cmd = 'BioSAK KEGG -t 4 -db_dir %s -seq_in %s -diamond' % (KEGG_db_dir, recipient_gene_same_ctg_gene_seqs)
# os.system(BioSAK_KEGG_cmd)
# os.system(BioSAK_COG_cmd)

# read in annotation results
gene_annot_C_dict = {}
gene_annot_D_dict = {}
for gene_annot in open(KEGG_annotation_results):
    gene_annot_split = gene_annot.strip().split('\t')
    gene_id = gene_annot_split[0]
    if len(gene_annot_split) == 1:
        gene_annot_C_dict[gene_id] = 'NA'
        gene_annot_D_dict[gene_id] = 'NA'
    else:
        if gene_id not in gene_annot_C_dict:
            gene_annot_C_dict[gene_id] = gene_annot_split[7]
            gene_annot_D_dict[gene_id] = gene_annot_split[8]
        else:
            if gene_annot_split[7] != gene_annot_C_dict[gene_id]:
                gene_annot_C_dict[gene_id] += '|%s' % gene_annot_split[7]
            if gene_annot_split[8] != gene_annot_D_dict[gene_id]:
                gene_annot_D_dict[gene_id] += '|%s' % gene_annot_split[8]

for recipient_gene in recipient_gene_list:

    co_ctg_gene_list = recipient_gene_same_ctg_genes_dict[recipient_gene]
    gene_num_in_ctg = len(co_ctg_gene_list)
    recipient_gene_pos = co_ctg_gene_list.index(recipient_gene)

    genes_on_the_left = []
    genes_on_the_right = []
    if (gene_num_in_ctg) > 1:
        if recipient_gene_pos > 0:
            genes_on_the_left = co_ctg_gene_list[0:recipient_gene_pos][::-1]
        if recipient_gene_pos < gene_num_in_ctg - 1:
            genes_on_the_right = co_ctg_gene_list[recipient_gene_pos + 1:]

    print(recipient_gene)

    # get the nearest HGT related functional gene in left side
    #print(genes_on_the_left)
    hgt_related_left_existence_D = []
    for left_side_gene in genes_on_the_left:
        left_side_gene_function_D = gene_annot_D_dict[left_side_gene]
        if check_list_elements_exist(key_word_list_D, left_side_gene_function_D) ==1:
            hgt_related_left_existence_D.append(left_side_gene)

    print(hgt_related_left_existence_D)


    # get the nearest HGT related functional gene in right side
    #print(genes_on_the_right)
    hgt_related_right_existence_D = []
    for right_side_gene in genes_on_the_right:
        right_side_gene_function_D = gene_annot_D_dict[right_side_gene]
        if check_list_elements_exist(key_word_list_D, right_side_gene_function_D) == 1:
            hgt_related_right_existence_D.append(right_side_gene)

    print(hgt_related_right_existence_D)
    print()






'''
module load python/3.7.3
source ~/mypython3env/bin/activate
module load diamond/0.9.24
cd /srv/scratch/z5039045
BioSAK KEGG -db_dir /srv/scratch/z5039045/DB/KEGG_2016-09-26 -t 12 -seq_in recipient_gene_same_ctg_genes.fasta -diamond
'''