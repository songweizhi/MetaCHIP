import os
__author__ = 'weizhisong'


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 30
walltime_needed = '05:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['']
output_folder = '/Users/songweizhi/Desktop/qsub_files/'

###########################################################################################


# Put '/' at the end of output_folder variant if not provided
if output_folder[-1] != '/':
    output_folder += '/'
else:
    pass

# Create output folder
os.system('mkdir ' + output_folder)

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae'
#line_8 = '\n\ncd $PBS_O_WORDIR\n'
line_8 = '\n\ncd /srv/scratch/z5095773/BP_genomes\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

################################################################

matches = open('/Users/songweizhi/Desktop/access.txt', 'r')
n = 0
for each in matches:
    each = each.strip()
    handle = open('/Users/songweizhi/Desktop/qsub_files/qsub_%s.sh' % each, 'w')
    handle.write(header)
    #handle.write(module_lines)
    handle.write('/srv/scratch/lanlab/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump --gzip --split-files %s' % each)
    handle.close()
    n += 1
print(n)
