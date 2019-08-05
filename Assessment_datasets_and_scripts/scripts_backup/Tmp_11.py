import os
import glob


wd = '/Users/songweizhi/Desktop/HgtSIM_wd/between_class'
coverage_cutoff = 75
min_taxa_num = 3


# class
all_vs_all_blast_results = 'selected_20_genomes_all_vs_all_ffn_class.tab'
donor_genome_marker = 'A'
recipient_genome_marker = 'G'
core_genes_id_file = 'core_genes_class_taxa3_cov80.txt'
core_genes_seq_file = 'core_genes_class_taxa3_cov80.fasta'
donor_ffn_folder = 'ffn_folder_alpha'
recipient_genome_folder = 'genome_folder_gamma'
recipient_genome_ext = 'fna'


# # genus
# all_vs_all_blast_results = 'selected_20_genomes_all_vs_all_ffn_genus.tab'
# donor_genome_marker = 'SB'
# recipient_genome_marker = 'SM'
# core_genes_id_file = 'core_genes_genus_taxa3_cov80.txt'
# core_genes_seq_file = 'core_genes_genus_taxa3_cov80.fasta'


os.chdir(wd)

# get
donor_ffn_file_re = '%s/*.ffn' % donor_ffn_folder
donor_ffn_file_list = sorted([os.path.basename(file_name) for file_name in glob.glob(donor_ffn_file_re)])
donor_genome_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in donor_ffn_file_list]

recipient_genome_re = '%s/*.%s' % (recipient_genome_folder, recipient_genome_ext)
recipient_genome_list = sorted([os.path.basename(file_name) for file_name in glob.glob(recipient_genome_re)])
recipient_genome_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in recipient_genome_list]



def get_genome_to_core_gene_dict(blast_results, donor_genomes, recipient_genomes, cov_cutoff, min_taxa_num):

    gene_ortholog_num_dict = {}
    for match in open(blast_results):
        match_split = match.strip().split('\t')
        query = match_split[0]
        subject = match_split[1]
        align_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        query_genome = '_'.join(query.split('_')[:-1])
        subject_genome = '_'.join(subject.split('_')[:-1])
        coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
        coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))
        if (query_genome in donor_genomes) and (query != subject):
            if (coverage_q >= cov_cutoff) and (coverage_s >= cov_cutoff):

                # initialize key value
                if query not in gene_ortholog_num_dict:
                    gene_ortholog_num_dict[query] = [set(), set()]

                if subject_genome in donor_genomes:
                    gene_ortholog_num_dict[query][0].add(subject_genome)

                if subject_genome in recipient_genomes:
                    gene_ortholog_num_dict[query][1].add(subject_genome)


    genome_qualified_gene_dict = {}
    for each_gene in gene_ortholog_num_dict:
        if (len(gene_ortholog_num_dict[each_gene][0]) >= min_taxa_num) and (len(gene_ortholog_num_dict[each_gene][1]) >= min_taxa_num):
            each_gene_genome = '_'.join(each_gene.split('_')[:-1])

            if each_gene_genome not in genome_qualified_gene_dict:
                genome_qualified_gene_dict[each_gene_genome] = {each_gene}
            else:
                genome_qualified_gene_dict[each_gene_genome].add(each_gene)

    return genome_qualified_gene_dict



donor_genome_core_gene_dict = get_genome_to_core_gene_dict(all_vs_all_blast_results, donor_genome_list_no_ext, recipient_genome_list_no_ext, coverage_cutoff, min_taxa_num)

print(len(genome_core_gene_dict))
