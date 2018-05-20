import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '02:59:00'
email = 'wythe1987@163.com'
modules_needed = []

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_gunzip'
wd_on_katana = '/srv/scratch/z5039045/hospital_effluents/raw_reads'

###########################################################################################

os.chdir(wd)

# create outputs folder
if not os.path.exists(outputs_folder):
    os.makedirs(outputs_folder)
else:
    shutil.rmtree(outputs_folder)
    os.makedirs(outputs_folder)

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'


#print(header)
#print(module_lines)

################################################################

for each in open('file_list_uniq.txt'):
    each = each.strip()
    output_handle = open('%s/qsub_gunzip_%s.sh' % (outputs_folder, each), 'w')
    reads_R1 = '%s_R1_001.fastq.gz' % each
    reads_R2 = '%s_R2_001.fastq.gz' % each
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('gunzip %s\n' % reads_R1)
    output_handle.write('gunzip %s\n' % reads_R2)
    output_handle.close()
