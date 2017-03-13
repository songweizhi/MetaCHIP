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

args = vars(parser.parse_args())
pwd_cfg_file = args['cfg']
config.read(pwd_cfg_file)

nodes_number = int(config['GENERAL']['qsub_nodes'])
ppn_number = int(config['GENERAL']['qsub_ppn'])
memory = int(config['STEP_4_BINNING']['qsub_memory_4'])
walltime_needed = config['STEP_4_BINNING']['qsub_walltime_4']
email = config['GENERAL']['qsub_email']
GemSIM_wd = config['GENERAL']['GemSIM_wd']
IDBA_UD_wd = config['GENERAL']['IDBA_UD_wd']
Mapping_wd = config['GENERAL']['Mapping_wd']
Binning_wd = config['GENERAL']['Binning_wd']
abundance_file = config['GENERAL']['abundance_file']
prefix = config['GENERAL']['prefix']
mink = int(config['STEP_2_ASSEMBLE']['idba_ud_mink'])
maxk = int(config['STEP_2_ASSEMBLE']['idba_ud_maxk'])
step = int(config['STEP_2_ASSEMBLE']['idba_ud_step'])
min_contig_length = int(config['STEP_3_MAPPING']['min_contig_length'])
pwd_jgi_summarize_bam_contig_depths = config['STEP_4_BINNING']['pwd_jgi_summarize_bam_contig_depths']
modules_metabat = config['STEP_4_BINNING']['qsub_modules_4_metabat']
modules_mycc = config['STEP_4_BINNING']['qsub_modules_4_mycc']

wd = os.getcwd()
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
modules_metabat_split = modules_metabat.split(',')
metabat_module_lines = ''
for module in modules_metabat_split:
    metabat_module_lines += 'module load ' + module + '\n'

modules_mycc_split = modules_mycc.split(',')
mycc_module_lines = ''
for module in modules_mycc_split:
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
qsub_mycc_file_handle.write('module unload python/3.5.2\n')
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
