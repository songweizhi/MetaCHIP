import os
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 3
memory = 30
walltime_needed = '02:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['idba/1.1.1']
genome_folder = 'input_genomes'
abundance_file = 'abundance.txt'
GemSIM_wd = '1_GemSIM'
IDBA_UD_wd = '2_IDBA_UD'
prefix = 'replicate'
mink = 20
maxk = 100
step = 20


#wd = os.getcwd()
wd = '/Users/songweizhi/Desktop/one_step'
pwd_genome_folder = '%s/%s' % (wd, genome_folder)
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)
pwd_IDBA_UD_wd = '%s/%s' % (wd, IDBA_UD_wd)

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
if os.path.isdir(pwd_IDBA_UD_wd):
    shutil.rmtree(pwd_IDBA_UD_wd)
    os.mkdir(pwd_IDBA_UD_wd)
else:
    os.mkdir(pwd_IDBA_UD_wd)

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

pwd_qsub_idba_ud_file = '%s/%s/qsub_idba_ud.sh' % (wd, 'qsub_files')
qsub_idba_ud_file_handle = open(pwd_qsub_idba_ud_file, 'w')

qsub_idba_ud_file_handle.write(header + module_lines)
qsub_idba_ud_file_handle.write('cat %s > %s\n' % (' '.join(R1_reads_files), pwd_combined_R1_fastq_file_name))
qsub_idba_ud_file_handle.write('cat %s > %s\n' % (' '.join(R2_reads_files), pwd_combined_R2_fastq_file_name))
qsub_idba_ud_file_handle.write('fq2fa --merge %s %s %s\n' % (pwd_combined_R1_fastq_file_name, pwd_combined_R2_fastq_file_name, pwd_combined_BP_fasta_file_name))
qsub_idba_ud_file_handle.write('idba_ud --pre_correction --num_threads %s --mink %s --maxk %s --step %s --read %s --out %s/combined_k%s-%s\n' % (ppn_number,
                                                                                                                                                 mink,
                                                                                                                                                 maxk,
                                                                                                                                                 step,
                                                                                                                                                 pwd_combined_BP_fasta_file_name,
                                                                                                                                                 pwd_IDBA_UD_wd,
                                                                                                                                                 mink,
                                                                                                                                                 maxk))

qsub_idba_ud_file_handle.close()
current_wd = os.getcwd()
os.chdir('%s/qsub_files' % wd)
os.system('qsub %s' % pwd_qsub_idba_ud_file)
os.chdir(current_wd)
