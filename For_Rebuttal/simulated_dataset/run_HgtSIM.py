import os
import glob
import random
import shutil
from Bio import SeqIO


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
        if (len(gene_ortholog_num_dict[each_gene][0]) >= min_taxa_num):
            each_gene_genome = '_'.join(each_gene.split('_')[:-1])

            if each_gene_genome not in genome_qualified_gene_dict:
                genome_qualified_gene_dict[each_gene_genome] = {each_gene}
            else:
                genome_qualified_gene_dict[each_gene_genome].add(each_gene)

    return genome_qualified_gene_dict


# common parameters
coverage_cutoff = 75
min_taxa_num = 3
transfer_per_genome = 10
mutation_levels = [0, 5, 10, 15, 20, 25, 30]
HgtSIM_script = '/Users/songweizhi/PycharmProjects/HgtSIM/HgtSIM.py'
flanking_seq = 'TAGATGAGTGATTAGTTAGTTA'


# # between class (A and G)
# wd = '/Users/songweizhi/Desktop/HgtSIM_wd/between_class_AG'
# all_vs_all_blast_results = 'selected_20_genomes_all_vs_all_ffn_class.tab'
# donor_ffn_folder = 'ffn_folder_alpha'
# recipient_genome_folder = 'genome_folder_gamma'
# recipient_genome_ext = 'fna'


# between class (A and B)
wd = '/Users/songweizhi/Desktop/HgtSIM_wd/between_class'
all_vs_all_blast_results = 'selected_20_genomes_all_vs_all_ffn_class.tab'
donor_ffn_folder = 'ffn_folder_alpha'
recipient_genome_folder = 'genome_folder_beta'
recipient_genome_ext = 'fna'


# # between genus
# wd = '/Users/songweizhi/Desktop/HgtSIM_wd/between_genus'
# all_vs_all_blast_results = 'selected_20_genomes_all_vs_all_ffn_genus.tab'
# donor_ffn_folder = 'ffn_folder_donor'
# recipient_genome_folder = 'genome_folder_SM'
# recipient_genome_ext = 'fna'


os.chdir(wd)


# Create qsub file folder
HgtSIM_input_files_folder = 'HgtSIM_input_files'
if os.path.isdir(HgtSIM_input_files_folder):
    shutil.rmtree(HgtSIM_input_files_folder)
    if os.path.isdir(HgtSIM_input_files_folder):
        shutil.rmtree(HgtSIM_input_files_folder)
        if os.path.isdir(HgtSIM_input_files_folder):
            shutil.rmtree(HgtSIM_input_files_folder)
os.system('mkdir %s' % HgtSIM_input_files_folder)

# get donor genome id list
donor_ffn_file_re = '%s/*.ffn' % donor_ffn_folder
donor_ffn_file_list = sorted([os.path.basename(file_name) for file_name in glob.glob(donor_ffn_file_re)])
donor_genome_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in donor_ffn_file_list]
print(donor_genome_list_no_ext)

# get recipient genome id list
recipient_genome_re = '%s/*.%s' % (recipient_genome_folder, recipient_genome_ext)
recipient_genome_list = sorted([os.path.basename(file_name) for file_name in glob.glob(recipient_genome_re)])
recipient_genome_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in recipient_genome_list]

# get donor genome core gene dict
donor_genome_core_gene_dict = get_genome_to_core_gene_dict(all_vs_all_blast_results, donor_genome_list_no_ext, recipient_genome_list_no_ext, coverage_cutoff, min_taxa_num)


for each in donor_genome_core_gene_dict:
    print(len(donor_genome_core_gene_dict[each]))

# run HgtSIM
bootstrap_num = 1
while bootstrap_num <= 10:
    print('Processing bootstrap%s' % bootstrap_num)

    # define output files
    pwd_seq_file = '%s/gene_seq_bootstrap_%s.fasta' % (HgtSIM_input_files_folder, bootstrap_num)
    distribution_file = '%s/gene_distribution_bootstrap_%s.txt' % (HgtSIM_input_files_folder, bootstrap_num)
    seq_file_handle = open(pwd_seq_file, 'w')
    selected_gene_dict = {}
    for donor_genome in donor_genome_list_no_ext:

        pwd_ffn_file = '%s/%s.ffn' % (donor_ffn_folder, donor_genome)
        gene_list_selected = random.sample(donor_genome_core_gene_dict[donor_genome], transfer_per_genome)
        selected_gene_dict[donor_genome] = gene_list_selected

        # export sequences
        for each_seq in SeqIO.parse(pwd_ffn_file, 'fasta'):
            if str(each_seq.id) in gene_list_selected:
                seq_file_handle.write('>%s\n' % str(each_seq.id))
                seq_file_handle.write('%s\n' % str(each_seq.seq))
    seq_file_handle.close()

    # get distribution file
    distribution_file_handle = open(distribution_file, 'w')
    n = 0
    for recipient_genome in recipient_genome_list_no_ext:
        for_write_list = [recipient_genome]
        for donor_genome in donor_genome_list_no_ext:
            for_write_list.append(selected_gene_dict[donor_genome][n])
        for_write = '%s\n' % (','.join(for_write_list))
        distribution_file_handle.write(for_write)
        n += 1
    distribution_file_handle.close()

    # prepare commands for HgtSIM
    output_prefix = 'bootstrap%s' % (bootstrap_num)
    for mutation_level in mutation_levels:
        HgtSIM_cmd = 'python3 %s -p %s -t %s -i %s -d %s -f %s -r 1-0-1-1 -x fna -lf %s -rf %s' % (HgtSIM_script, output_prefix, pwd_seq_file, mutation_level, distribution_file, recipient_genome_folder, flanking_seq, flanking_seq)
        os.system(HgtSIM_cmd)

    # Create current bootstrap folder
    if os.path.isdir(output_prefix):
        shutil.rmtree(output_prefix)
        if os.path.isdir(output_prefix):
            shutil.rmtree(output_prefix)
            if os.path.isdir(output_prefix):
                shutil.rmtree(output_prefix)

    os.system('mkdir %s' % output_prefix)
    os.system('mv %s_* %s/' % (output_prefix, output_prefix))

    bootstrap_num += 1

