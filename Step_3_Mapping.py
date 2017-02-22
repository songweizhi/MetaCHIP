import os
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 3
memory = 30
walltime_needed = '02:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['bowtie/2.2.9', 'samtools/1.2']
genome_folder = 'input_genomes'
abundance_file = 'abundance.txt'
GemSIM_wd = '1_GemSIM'
IDBA_UD_wd = '2_IDBA_UD'
Mapping_wd = '3_Mapping'
prefix = 'replicate'
mink = 20
maxk = 100
step = 20
pwd_select_contig_pl = '/srv/scratch/z5039045/Scripts/select_contig.pl'
pwd_get_fasta_stats_pl = '/srv/scratch/z5039045/Scripts/get_fasta_stats.pl'
min_contig_length = 2500

wd = os.getcwd()
#wd = '/Users/songweizhi/Desktop/one_step'
pwd_genome_folder = '%s/%s' % (wd, genome_folder)
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)
pwd_IDBA_UD_wd = '%s/%s' % (wd, IDBA_UD_wd)
pwd_Mapping_wd = '%s/%s' % (wd, Mapping_wd)

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
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

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
if os.path.isdir(pwd_Mapping_wd):
    shutil.rmtree(pwd_Mapping_wd)
    os.mkdir(pwd_Mapping_wd)
else:
    os.mkdir(pwd_Mapping_wd)

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

os.system('perl %s -m %s %s %s' % (pwd_select_contig_pl, min_contig_length, pwd_scaffold, pwd_idba_ud_output_scaffold_filtered))
os.system('perl %s -T %s > %s.txt' % (pwd_get_fasta_stats_pl, pwd_idba_ud_output_scaffold_filtered, pwd_idba_ud_output_scaffold_filtered))
os.system('module load %s' % 'bowtie/2.2.9')
os.system('bowtie2-build -f %s %s' % (pwd_idba_ud_output_scaffold_filtered, pwd_idba_ud_output_scaffold_filtered_no_extension))

n = 1
while n <= number_of_replicates:
    pwd_qsub_mapping_file = '%s/%s/qsub_mapping_%s_%s.sh' % (wd, 'qsub_files', prefix, n)
    qsub_mapping_file_handle = open(pwd_qsub_mapping_file, 'w')
    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_GemSIM_wd, n)
    fastq_file_R1_Q30_P = '%s_%s_fir_Q30_P.fastq' % (prefix, n)
    fastq_file_R2_Q30_P = '%s_%s_sec_Q30_P.fastq' % (prefix, n)
    pwd_fastq_file_R1_Q30_P = '%s/%s' % (pwd_current_replicate_wd, fastq_file_R1_Q30_P)
    pwd_fastq_file_R2_Q30_P = '%s/%s' % (pwd_current_replicate_wd, fastq_file_R2_Q30_P)
    pwd_sam_file = '%s/%s_%s.sam' % (pwd_Mapping_wd, prefix, n)
    pwd_bam_file = '%s/%s_%s.bam' % (pwd_Mapping_wd, prefix, n)
    pwd_bam_file_sorted = '%s/%s_%s_sorted' % (pwd_Mapping_wd, prefix, n)
    qsub_mapping_file_handle.write(header + module_lines)
    qsub_mapping_file_handle.write('bowtie2 -x %s -1 %s -2 %s -S %s -p 6 -q\n' % (pwd_idba_ud_output_scaffold_filtered_no_extension, pwd_fastq_file_R1_Q30_P, pwd_fastq_file_R2_Q30_P, pwd_sam_file))
    qsub_mapping_file_handle.write('samtools view -bS %s -o %s\n' % (pwd_sam_file, pwd_bam_file))
    qsub_mapping_file_handle.write('samtools sort %s %s\n' % (pwd_bam_file, pwd_bam_file_sorted))
    qsub_mapping_file_handle.write('samtools index %s.bam\n' % pwd_bam_file_sorted)
    qsub_mapping_file_handle.close()

    current_wd = os.getcwd()
    os.chdir('%s/qsub_files' % wd)
    os.system('qsub %s' % pwd_qsub_mapping_file)
    os.chdir(current_wd)

    n += 1
