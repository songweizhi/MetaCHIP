import os
import shutil
import argparse
import configparser

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 3
memory = 30
walltime_needed = '02:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_metabat = ['metabat/0.32.4']
modules_mycc = ['intel/16.0.1.150', 'perl/5.20.1', 'python/2.7.10', 'cd-hit/4.6.4', 'prodigal/2.6.3', 'parallel/20160222', 'hmmer/3.1b2', 'barrnap/0.7', 'mycc/20150710']

genome_folder = 'input_genomes'
abundance_file = 'abundance.txt'
GemSIM_wd = '1_GemSIM'
IDBA_UD_wd = '2_IDBA_UD'
Mapping_wd = '3_Mapping'
Binning_wd = '4_Binning'
prefix = 'replicate'
mink = 20
maxk = 100
step = 20
min_contig_length = 2500
pwd_jgi_summarize_bam_contig_depths = '/share/apps/metabat/0.32.4/jgi_summarize_bam_contig_depths'

wd = os.getcwd()
#wd = '/Users/songweizhi/Desktop/one_step'
pwd_genome_folder = '%s/%s' % (wd, genome_folder)
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)
pwd_IDBA_UD_wd = '%s/%s' % (wd, IDBA_UD_wd)
pwd_Mapping_wd = '%s/%s' % (wd, Mapping_wd)
pwd_Binning_wd = '%s/%s' % (wd, Binning_wd)
pwd_qsub_file_folder = '%s/qsub_files' % wd

###################################### Prepare qsub file header ######################################

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae\n\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
metabat_module_lines = ''
for module in modules_metabat:
    metabat_module_lines += 'module load ' + module + '\n'

mycc_module_lines = ''
for module in modules_mycc:
    mycc_module_lines += 'module load ' + module + '\n'

######################################## Prepare input files ########################################

# get combined abundance list:
abundance_combined = open(pwd_abundance_file)
combined_abundance_list = []
for each in abundance_combined:
    each_split = each.strip().split('\t')
    replicate_number = len(each_split) - 1
    combined_abundance_list.append(each_split)
number_of_replicates = len(combined_abundance_list[0]) - 1


# prepare GemSIM working directory
if os.path.isdir(pwd_Binning_wd):
    shutil.rmtree(pwd_Binning_wd)
    os.mkdir(pwd_Binning_wd)
else:
    os.mkdir(pwd_Binning_wd)

n = 1
R1_reads_files = []
R2_reads_files = []
while n <= number_of_replicates:
    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_GemSIM_wd, n)
    fastq_file_R1_Q30_P = '%s_%s_fir_Q30_P.fastq' % (prefix, n)
    fastq_file_R2_Q30_P = '%s_%s_sec_Q30_P.fastq' % (prefix, n)
    cmd_fastqc = 'fastqc %s %s\n' % (fastq_file_R1_Q30_P, fastq_file_R2_Q30_P)
    pwd_fastq_file_R1_Q30_P = '%s/%s' % (pwd_current_replicate_wd, fastq_file_R1_Q30_P)
    pwd_fastq_file_R2_Q30_P = '%s/%s' % (pwd_current_replicate_wd, fastq_file_R2_Q30_P)
    R1_reads_files.append(pwd_fastq_file_R1_Q30_P)
    R2_reads_files.append(pwd_fastq_file_R2_Q30_P)
    n += 1

combined_R1_fastq_file_name = 'combined_R1.fastq'
combined_R2_fastq_file_name = 'combined_R2.fastq'
combined_BP_fasta_file_name = 'combined.fasta'
pwd_combined_R1_fastq_file_name = '%s/%s' % (pwd_IDBA_UD_wd, combined_R1_fastq_file_name)
pwd_combined_R2_fastq_file_name = '%s/%s' % (pwd_IDBA_UD_wd, combined_R2_fastq_file_name)
pwd_combined_BP_fasta_file_name = '%s/%s' % (pwd_IDBA_UD_wd, combined_BP_fasta_file_name)


idba_ud_output_scaffold = 'scaffold.fa'
idba_ud_output_scaffold_filtered = 'scaffold_k%s-%s_lt%s.fa' % (mink, maxk, min_contig_length)
idba_ud_output_scaffold_filtered_no_extension = 'scaffold_k%s-%s_lt%s' % (mink, maxk, min_contig_length)
pwd_scaffold = '%s/combined_k%s-%s/%s' % (pwd_IDBA_UD_wd, mink, maxk, idba_ud_output_scaffold)
pwd_idba_ud_output_scaffold_filtered = '%s/%s' % (pwd_Mapping_wd, idba_ud_output_scaffold_filtered)
pwd_idba_ud_output_scaffold_filtered_no_extension = '%s/%s' % (pwd_Mapping_wd, idba_ud_output_scaffold_filtered_no_extension)


n = 1
sorted_bam_files = []
while n <= number_of_replicates:
    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_GemSIM_wd, n)
    fastq_file_R1_Q30_P = '%s_%s_fir_Q30_P.fastq' % (prefix, n)
    fastq_file_R2_Q30_P = '%s_%s_sec_Q30_P.fastq' % (prefix, n)
    pwd_fastq_file_R1_Q30_P = '%s/%s' % (pwd_current_replicate_wd, fastq_file_R1_Q30_P)
    pwd_fastq_file_R2_Q30_P = '%s/%s' % (pwd_current_replicate_wd, fastq_file_R2_Q30_P)
    pwd_bam_file = '%s/%s_%s.bam' % (pwd_Mapping_wd, prefix, n)
    pwd_sorted_bam_file = '%s/%s_%s_sorted.bam' % (pwd_Mapping_wd, prefix, n)
    sorted_bam_files.append(pwd_sorted_bam_file)
    n += 1


pwd_depth_file = '%s/%s_depth.txt' % (pwd_Binning_wd, idba_ud_output_scaffold_filtered_no_extension)
pwd_paired_contig_file = '%s/%s_paired.txt' % (pwd_Binning_wd, idba_ud_output_scaffold_filtered_no_extension)
os.system('%s --outputDepth %s --pairedContigs %s %s' % (pwd_jgi_summarize_bam_contig_depths, pwd_depth_file, pwd_paired_contig_file, ' '.join(sorted_bam_files)))
metabat_bin_folder = '%s/%s_MetaBAT' % (pwd_Binning_wd, idba_ud_output_scaffold_filtered_no_extension)


pwd_qsub_metabat_file = '%s/qsub_MetaBAT.sh' % (pwd_qsub_file_folder)
pwd_qsub_mycc_file = '%s/qsub_MyCC.sh' % (pwd_qsub_file_folder)
qsub_metabat_file_handle = open(pwd_qsub_metabat_file, 'w')
qsub_mycc_file_handle = open(pwd_qsub_mycc_file, 'w')


qsub_metabat_file_handle.write(header + metabat_module_lines)
qsub_metabat_file_handle.write('cd %s\n' % pwd_Binning_wd)
qsub_metabat_file_handle.write('mkdir %s\n' % metabat_bin_folder)
qsub_metabat_file_handle.write('metabat -i %s -a %s -p %s -o %s/%s_MetaBAT/%s' % (pwd_idba_ud_output_scaffold_filtered, pwd_depth_file, pwd_paired_contig_file, pwd_Binning_wd, idba_ud_output_scaffold_filtered_no_extension, idba_ud_output_scaffold_filtered_no_extension))
qsub_metabat_file_handle.close()
qsub_mycc_file_handle.write(header)
qsub_mycc_file_handle.write('module unload intel/11.1.080\n')
qsub_mycc_file_handle.write(mycc_module_lines)
qsub_mycc_file_handle.write('\ncd %s\n' % pwd_Binning_wd)
qsub_mycc_file_handle.write('MyCC.py %s -a %s 56mer' % (pwd_idba_ud_output_scaffold_filtered, pwd_depth_file))
qsub_mycc_file_handle.close()


current_wd = os.getcwd()
os.chdir('%s/qsub_files' % wd)
os.system('qsub %s' % pwd_qsub_metabat_file)
os.system('qsub -V %s' % pwd_qsub_mycc_file)
os.chdir(current_wd)
n += 1
