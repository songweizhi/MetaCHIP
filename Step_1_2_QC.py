import os

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 3
memory = 30
walltime_needed = '00:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['java/8u91', 'fastqc/0.10.1']
genome_folder = 'input_genomes'
abundance_file = 'abundance.txt'
GemSIM_wd = '1_GemSIM'
prefix = 'replicate'

pwd_trimmomatic_executable = '/share/apps/trimmomatic/0.33/trimmomatic-0.33.jar'
trimmomatic_parameters = 'ILLUMINACLIP:/share/apps/trimmomatic/0.33/adapters/TruSeq3-PE-2.fa:2:30:10:6:true LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:50'

wd = os.getcwd()
#wd = '/Users/songweizhi/Desktop/one_step'
pwd_genome_folder = '%s/%s' % (wd, genome_folder)
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)

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

# get number of replicates:
abundance_combined = open(pwd_abundance_file)
combined_abundance_list = []
for each in abundance_combined:
    each_split = each.strip().split('\t')
    combined_abundance_list.append(each_split)
number_of_replicates = len(combined_abundance_list[0]) - 1

n = 1
while n <= number_of_replicates:
    pwd_qsub_qc_file = '%s/%s/%s/qsub_QC_%s_%s.sh' % (wd, GemSIM_wd, 'qsub_files', prefix, n)
    qsub_qc_handle = open(pwd_qsub_qc_file, 'w')
    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_GemSIM_wd, n)
    fastq_file_R1 = '%s_%s_fir.fastq' % (prefix, n)
    fastq_file_R2 = '%s_%s_sec.fastq' % (prefix, n)
    fastq_file_R1_Q30_P = '%s_%s_fir_Q30_P.fastq' % (prefix, n)
    fastq_file_R1_Q30_UP = '%s_%s_fir_Q30_UP.fastq' % (prefix, n)
    fastq_file_R2_Q30_P = '%s_%s_sec_Q30_P.fastq' % (prefix, n)
    fastq_file_R2_Q30_UP = '%s_%s_sec_Q30_UP.fastq' % (prefix, n)
    cmd_trimmomatic = 'java -jar %s PE %s %s %s %s %s %s %s\n' % (pwd_trimmomatic_executable,
                                                                fastq_file_R1,
                                                                fastq_file_R2,
                                                                fastq_file_R1_Q30_P,
                                                                fastq_file_R1_Q30_UP,
                                                                fastq_file_R2_Q30_P,
                                                                fastq_file_R2_Q30_UP,
                                                                trimmomatic_parameters)

    cmd_fastqc = 'fastqc %s %s\n' % (fastq_file_R1_Q30_P, fastq_file_R2_Q30_P)
    qsub_qc_handle.write(header + module_lines)
    qsub_qc_handle.write('cd %s\n' % pwd_current_replicate_wd)
    qsub_qc_handle.write(cmd_trimmomatic)
    qsub_qc_handle.write(cmd_fastqc)
    qsub_qc_handle.close()

    current_wd = os.getcwd()
    os.chdir('%s/qsub_files' % pwd_GemSIM_wd)
    os.system('qsub %s' % pwd_qsub_qc_file)
    os.chdir(current_wd)
    n += 1
