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
memory = int(config['STEP_1_1_GemSIM']['qsub_memory_1_1'])
walltime_needed = config['STEP_1_1_GemSIM']['qsub_walltime_1_1']
email = config['GENERAL']['qsub_email']
modules_needed = config['STEP_1_1_GemSIM']['qsub_modules_1_1']
genome_folder = config['GENERAL']['genome_folder']
abundance_file = config['GENERAL']['abundance_file']
GemSIM_wd = config['GENERAL']['GemSIM_wd']
prefix = config['GENERAL']['prefix']
pwd_GemReads_executable = config['STEP_1_1_GemSIM']['pwd_GemReads_executable']
pwd_error_models = config['STEP_1_1_GemSIM']['pwd_error_models']
GemReads_parameters = config['STEP_1_1_GemSIM']['GemReads_parameters']

wd = os.getcwd()
pwd_genome_folder = '%s/%s' % (wd, genome_folder)
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

# get combined abundance list:
abundance_combined = open(pwd_abundance_file)
combined_abundance_list = []
for each in abundance_combined:
    each_split = each.strip().split('\t')
    replicate_number = len(each_split) - 1
    combined_abundance_list.append(each_split)
number_of_genome = len(combined_abundance_list)
number_of_replicates = len(combined_abundance_list[0]) - 1


# prepare GemSIM working directory
if os.path.isdir(pwd_GemSIM_wd):
    shutil.rmtree(pwd_GemSIM_wd)
    os.mkdir(pwd_GemSIM_wd)
else:
    os.mkdir(pwd_GemSIM_wd)


if os.path.isdir(pwd_qsub_file_folder):
    shutil.rmtree(pwd_qsub_file_folder)
    os.mkdir(pwd_qsub_file_folder)
else:
    os.mkdir(pwd_qsub_file_folder)


n = 1
while n <= number_of_replicates:
    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_GemSIM_wd, n)
    os.mkdir(pwd_current_replicate_wd)
    # prepare current abundance file
    pwd_current_abundance_file = '%s/abundance_%s.txt' % (pwd_current_replicate_wd, n)
    current_abundance_handle = open(pwd_current_abundance_file, 'w')
    for each in combined_abundance_list:
        current_abundance_handle.write('%s\t%s\n' % (each[0], each[n]))

    # prepare current qsub file
    pwd_current_qsub_file = '%s/qsub_files/qsub_GemSIM_abundance_%s.sh' % (wd, n)
    current_qsub_handle = open(pwd_current_qsub_file, 'w')
    cmd = 'python %s -R ../../%s -a %s -m %s -o %s_%s %s\n' % (pwd_GemReads_executable,
                                                               genome_folder,
                                                               pwd_current_abundance_file,
                                                               pwd_error_models,
                                                               prefix,
                                                               n,
                                                               GemReads_parameters)

    current_qsub_handle.write(header + module_lines)
    current_qsub_handle.write('cd %s\n' % (pwd_current_replicate_wd))
    current_qsub_handle.write(cmd)
    current_abundance_handle.close()
    current_qsub_handle.close()

    current_wd = os.getcwd()
    os.chdir(pwd_qsub_file_folder)
    os.system('qsub %s' % pwd_current_qsub_file)
    os.chdir(current_wd)
    n += 1
