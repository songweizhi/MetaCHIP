import os
import shutil
import argparse
import configparser

###################################### CONFIGURATION ######################################

parser = argparse.ArgumentParser()
config = configparser.ConfigParser()

parser.add_argument('-cfg',
                    required=True,
                    help='path to configuration file')

parser.add_argument('-qsub',
                    action="store_true",
                    help='submit job script to the queue')

args = vars(parser.parse_args())
pwd_cfg_file = args['cfg']
config.read(pwd_cfg_file)

nodes_number = int(config['GENERAL']['qsub_nodes'])
ppn_number = int(config['STEP_2_ASSEMBLE_metaSPAdes']['qsub_ppn_2'])
qsub_folder = config['GENERAL']['qsub_folder']
memory = int(config['STEP_2_ASSEMBLE_metaSPAdes']['qsub_memory_2'])
walltime_needed = config['STEP_2_ASSEMBLE_metaSPAdes']['qsub_walltime_2']
email = config['GENERAL']['qsub_email']
modules_needed = config['STEP_2_ASSEMBLE_metaSPAdes']['qsub_modules_2']
abundance_file = config['GENERAL']['abundance_file']
GemSIM_wd = config['GENERAL']['GemSIM_wd']
metaSPAdes_wd = config['GENERAL']['metaSPAdes_wd']
prefix = config['GENERAL']['prefix']
kmer_range = config['STEP_2_ASSEMBLE_metaSPAdes']['kmer_range']

wd = os.getcwd()
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)
pwd_metaSPAdes_wd = '%s/%s' % (wd, metaSPAdes_wd)

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

# get combined abundance list:
abundance_combined = open(pwd_abundance_file)
combined_abundance_list = []
for each in abundance_combined:
    each_split = each.strip().split('\t')
    replicate_number = len(each_split) - 1
    combined_abundance_list.append(each_split)
number_of_replicates = len(combined_abundance_list[0]) - 1


# prepare metaSPAdes working directory
if os.path.isdir(pwd_metaSPAdes_wd):
    shutil.rmtree(pwd_metaSPAdes_wd)
    os.mkdir(pwd_metaSPAdes_wd)
else:
    os.mkdir(pwd_metaSPAdes_wd)

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
pwd_combined_R1_fastq_file_name = '%s/%s' % (pwd_metaSPAdes_wd, combined_R1_fastq_file_name)
pwd_combined_R2_fastq_file_name = '%s/%s' % (pwd_metaSPAdes_wd, combined_R2_fastq_file_name)

pwd_qsub_metaSPAdes_file = '%s/%s/qsub_metaSPAdes.sh' % (wd, qsub_folder)
qsub_metaSPAdes_file_handle = open(pwd_qsub_metaSPAdes_file, 'w')

qsub_metaSPAdes_file_handle.write(header + module_lines)
qsub_metaSPAdes_file_handle.write('cat %s > %s\n' % (' '.join(R1_reads_files), pwd_combined_R1_fastq_file_name))
qsub_metaSPAdes_file_handle.write('cat %s > %s\n' % (' '.join(R2_reads_files), pwd_combined_R2_fastq_file_name))
qsub_metaSPAdes_file_handle.write('spades.py --meta -t %s -1 %s -2 %s -o %s/combined_k21-127 -k %s\n' % (ppn_number, pwd_combined_R1_fastq_file_name, pwd_combined_R2_fastq_file_name, pwd_metaSPAdes_wd, kmer_range))

qsub_metaSPAdes_file_handle.close()

current_wd = os.getcwd()
os.chdir('%s/%s' % (wd, qsub_folder))
os.system('qsub %s' % pwd_qsub_metaSPAdes_file)
os.chdir(current_wd)
