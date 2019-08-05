import os
import shutil


##################################################### CONFIGURATION ####################################################

# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 12
memory = 60
walltime_needed = '11:59:00'
email = '244289990@qq.com'
modules_needed = ['python/2.7.12', 'blast+/2.6.0']
qsub_file_folder = '/Users/songweizhi/Desktop/qsub_files_BM'
#wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/between_genus/MetaCHIP_wd'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/between_class/MetaCHIP_wd'
mutation_levels = [0, 5, 10, 15, 20, 25, 30]


################################################## generate_qsub_file ##################################################


# Create qsub file folder
if os.path.isdir(qsub_file_folder):
    shutil.rmtree(qsub_file_folder)
    if os.path.isdir(qsub_file_folder):
        shutil.rmtree(qsub_file_folder)
        if os.path.isdir(qsub_file_folder):
            shutil.rmtree(qsub_file_folder)
os.system('mkdir %s' % qsub_file_folder)


# Prepare header
line_1 = '#!/bin/bash'
line_2 = '#PBS -l nodes=%s:ppn=%s' % (str(nodes_number), str(ppn_number))
line_3 = '#PBS -l vmem=%sgb' % str(memory)
line_4 = '#PBS -l walltime=%s' % walltime_needed
line_5 = '#PBS -j oe'
line_6 = '#PBS -M %s' % email
line_7 = '#PBS -m ae'
header = '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % (line_1, line_2, line_3, line_4, line_5, line_6, line_7)


# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load %s\n' % module


# write to qsub files
bootstrap_num = 1
while bootstrap_num <= 10:
    for mutation_level in mutation_levels:
        current_qsub_file = '%s/qsub_BM_bootstrap%s_m%s.sh' % (qsub_file_folder, bootstrap_num, mutation_level)
        handle = open(current_qsub_file, 'w')
        handle.write('\n' + header + '\n')
        handle.write(module_lines)
        handle.write('\ncd %s/bootstrap%s/bootstrap%s_m%s\n' % (wd_on_katana, bootstrap_num, bootstrap_num, mutation_level) )
        #BM_cmd = 'python /srv/scratch/z5039045/Softwares/MetaCHIP/Best-match.py -p bootstrap%s_m%s -g bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_grouping_g2.txt -threads 12' % (bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        BM_cmd = 'python /srv/scratch/z5039045/Softwares/MetaCHIP/Best-match.py -p bootstrap%s_m%s -g bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_grouping_c2.txt -threads 12' % (bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        handle.write(BM_cmd)
        handle.close()
    bootstrap_num += 1

