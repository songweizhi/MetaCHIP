import os
import shutil


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 12
memory = 120
walltime_needed = '23:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['python/2.7.9']


wd = os.getcwd()
pwd_GemReads_executable = '/srv/scratch/z5039045/Softwares/GemSIM_v1.6/GemReads.py'
pwd_genome_folder = '/srv/scratch/z5039045/MetaCHIP_wd/Simulated_datasets/30_genomes'
pwd_abundance_file = '/srv/scratch/z5039045/MetaCHIP_wd/Simulated_datasets/abundance.txt'
pwd_error_models = '/srv/scratch/z5039045/Softwares/GemSIM_v1.6/models/ill100v5_p.gzip'
GemReads_parameters = '-q 33 -u d -s 30 -n 10000000 -l d -p'
prefix = 'simulation_original'
#wd = '/Users/songweizhi/Desktop/test'

#pwd_abundance_file =    '%s/abundance.txt'  %   wd
#pwd_genome_folder =     '%s/30_genomes'     %   wd



###################################### Prepare qsub file header ######################################

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae'
line_8 = '\n\ncd $PBS_O_WORDIR\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'


# get separated abundance file:
abundance_combined = open(pwd_abundance_file)
n = 1
replicate_number = 0
abundance_list = []
for each in abundance_combined:
    each_split = each.strip().split('\t')
    replicate_number = len(each_split) - 1
    while n <= replicate_number:
        current_abundance = [each_split[0], each_split[n]]
        abundance_list.append(current_abundance)
        n += 1
    n = 1

abundance_file_path, abundance_file_name = os.path.split(pwd_abundance_file)
abundance_file_name_no_ext = '.'.join(abundance_file_name.split('.')[:-1])
abundance_file_name_ext = abundance_file_name.split('.')[-1]
characters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']

# prepare folders
if os.path.isdir('%s/GemSIM' % abundance_file_path):
    shutil.rmtree('%s/GemSIM' % abundance_file_path)
    os.mkdir('%s/GemSIM' % abundance_file_path)
else:
    os.mkdir('%s/GemSIM' % abundance_file_path)

m = 0
while m <= replicate_number - 1:
    os.mkdir('%s/GemSIM/%s_%s' % (abundance_file_path, prefix, characters[m]))
    pwd_abundance_file_sep = '%s/GemSIM/%s_%s/%s_%s.%s' % (abundance_file_path, prefix, characters[m], abundance_file_name_no_ext, characters[m], abundance_file_name_ext)
    abundance_sep = open(pwd_abundance_file_sep, 'w')
    current = m
    while current <= len(abundance_list) - 1:
        abundance_sep.write('%s\n' % '\t'.join(abundance_list[current]))
        current += replicate_number
    abundance_sep.close()
    command = 'python %s -R %s -a %s -m %s -o %s_%s %s' % (pwd_GemReads_executable,
                                                           pwd_genome_folder,
                                                           pwd_abundance_file_sep,
                                                           pwd_error_models,
                                                           prefix,
                                                           characters[m],
                                                           GemReads_parameters)

    pwd_qsub_file = '%s/GemSIM/%s_%s/qsub_%s_%s.sh' % (abundance_file_path, prefix, characters[m], abundance_file_name_no_ext, characters[m])
    qsub_file_handle = open(pwd_qsub_file, 'w')

    qsub_file_handle.write(header + module_lines)
    qsub_file_handle.write('%s\n' % command)

    qsub_file_handle.close()


    # cd to current folder
    #os.chdir('%s/GemSIM/%s_%s' % (abundance_file_path, prefix, characters[m]))

    # submit job
    #os.system('qsub %s' % pwd_qsub_file)
    print('qsub %s' % pwd_qsub_file)

    # cd back to wd
    os.chdir(wd)

    m += 1





