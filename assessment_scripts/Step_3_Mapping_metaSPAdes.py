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
qsub_folder = config['GENERAL']['qsub_folder']
memory = int(config['STEP_3_MAPPING_metaSPAdes']['qsub_memory_3'])
walltime_needed = config['STEP_3_MAPPING_metaSPAdes']['qsub_walltime_3']
email = config['GENERAL']['qsub_email']
modules_needed = config['STEP_3_MAPPING_metaSPAdes']['qsub_modules_3']
abundance_file = config['GENERAL']['abundance_file']
GemSIM_wd = config['GENERAL']['GemSIM_wd']
metaSPAdes_wd = config['GENERAL']['metaSPAdes_wd']
Mapping_wd = config['GENERAL']['Mapping_wd']
prefix = config['GENERAL']['prefix']
pwd_select_contig_pl = config['STEP_3_MAPPING_metaSPAdes']['pwd_select_contig_pl']
pwd_get_fasta_stats_pl = config['STEP_3_MAPPING_metaSPAdes']['pwd_get_fasta_stats_pl']
pwd_bowtie2_build = config['STEP_3_MAPPING_metaSPAdes']['pwd_bowtie2_build']
min_contig_length = int(config['STEP_3_MAPPING_metaSPAdes']['min_contig_length'])

wd = os.getcwd()
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_GemSIM_wd = '%s/%s' % (wd, GemSIM_wd)
pwd_metaSPAdes_wd = '%s/%s' % (wd, metaSPAdes_wd)
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
pwd_combined_R1_fastq_file_name = '%s/%s' % (pwd_metaSPAdes_wd, combined_R1_fastq_file_name)
pwd_combined_R2_fastq_file_name = '%s/%s' % (pwd_metaSPAdes_wd, combined_R2_fastq_file_name)

metaSPAdes_output_scaffold = 'scaffolds_renamed.fasta'
metaSPAdes_output_scaffold_filtered = 'scaffolds_k21-127_lt%s.fa' % (min_contig_length)
metaSPAdes_output_scaffold_filtered_no_extension = 'scaffolds_k21-127_lt%s' % (min_contig_length)
pwd_scaffold = '%s/combined_k21-127/%s' % (pwd_metaSPAdes_wd, metaSPAdes_output_scaffold)
pwd_metaSPAdes_output_scaffold_filtered = '%s/%s' % (pwd_Mapping_wd, metaSPAdes_output_scaffold_filtered)
pwd_metaSPAdes_output_scaffold_filtered_no_extension = '%s/%s' % (pwd_Mapping_wd, metaSPAdes_output_scaffold_filtered_no_extension)

os.system('perl %s -m %s %s %s' % (pwd_select_contig_pl, min_contig_length, pwd_scaffold, pwd_metaSPAdes_output_scaffold_filtered))
os.system('perl %s -T %s > %s.txt' % (pwd_get_fasta_stats_pl, pwd_metaSPAdes_output_scaffold_filtered, pwd_metaSPAdes_output_scaffold_filtered))
os.system('%s -f %s %s' % (pwd_bowtie2_build, pwd_metaSPAdes_output_scaffold_filtered, pwd_metaSPAdes_output_scaffold_filtered_no_extension))

n = 1
while n <= number_of_replicates:
    pwd_qsub_mapping_file = '%s/%s/qsub_mapping_%s_%s.sh' % (wd, qsub_folder, prefix, n)
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
    qsub_mapping_file_handle.write('bowtie2 -x %s -1 %s -2 %s -S %s -p 1 -q\n' % (pwd_metaSPAdes_output_scaffold_filtered_no_extension, pwd_fastq_file_R1_Q30_P, pwd_fastq_file_R2_Q30_P, pwd_sam_file))
    qsub_mapping_file_handle.write('samtools view -bS %s -o %s\n' % (pwd_sam_file, pwd_bam_file))
    qsub_mapping_file_handle.write('samtools sort %s %s\n' % (pwd_bam_file, pwd_bam_file_sorted))
    qsub_mapping_file_handle.write('samtools index %s.bam\n' % pwd_bam_file_sorted)
    qsub_mapping_file_handle.write('rm %s\n' % pwd_sam_file)
    qsub_mapping_file_handle.write('rm %s\n' % pwd_bam_file)

    qsub_mapping_file_handle.close()
    current_wd = os.getcwd()
    os.chdir('%s/%s' % (wd, qsub_folder))
    os.system('qsub %s' % pwd_qsub_mapping_file)
    os.chdir(current_wd)
    n += 1
