import os
import argparse
import configparser

###################################### CONFIGURATION ######################################

parser = argparse.ArgumentParser()
config = configparser.ConfigParser()

parser.add_argument('-cfg',
                    required=True,
                    help='path to configuration file')

args = vars(parser.parse_args())
pwd_cfg_file = args['cfg']
config.read(pwd_cfg_file)

nodes_number = int(config['GENERAL']['qsub_nodes'])
ppn_number = int(config['GENERAL']['qsub_ppn'])
memory = int(config['STEP_1_2_QC']['qsub_memory_1_2'])
walltime_needed = config['STEP_1_2_QC']['qsub_walltime_1_2']
email = config['GENERAL']['qsub_email']
modules_needed = config['STEP_1_2_QC']['qsub_modules_1_2']
abundance_file = config['GENERAL']['abundance_file']
GemSIM_wd = config['GENERAL']['GemSIM_wd']
prefix = config['GENERAL']['prefix']
pwd_trimmomatic_executable = config['STEP_1_2_QC']['pwd_trimmomatic_executable']
trimmomatic_parameters = config['STEP_1_2_QC']['trimmomatic_parameters']

wd = os.getcwd()
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)
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
modules_needed_split = modules_needed.split(',')
module_lines = ''
for module in modules_needed_split:
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
    pwd_qsub_qc_file = '%s/qsub_QC_%s_%s.sh' % (pwd_qsub_file_folder, prefix, n)
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
    os.chdir(pwd_qsub_file_folder)
    os.system('qsub %s' % pwd_qsub_qc_file)
    os.chdir(current_wd)
    n += 1
