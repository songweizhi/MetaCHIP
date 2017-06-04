import os
__author__ = 'weizhisong'


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 6
memory = 126
walltime_needed = '23:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['idba/1.1.1_512']
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
#line_8 = '\n\ncd /srv/scratch/z5095773/BP_genomes\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

################################################################

for each in open('/Users/songweizhi/Desktop/test.txt'):
    each_split = each.strip().split('/')
    input_fasta = '/srv/scratch/z5039045/HgtSIM/reads_length/%s/combined.fasta' % (each.strip())
    handle = open('/Users/songweizhi/Desktop/qsub_files/qsub_%s_%s.sh' % (each_split[0], each_split[1]), 'w')
    handle.write('\n' + header + '\n')
    handle.write(module_lines)
    handle.write('\ncd /srv/scratch/z5039045/HgtSIM/reads_length/%s\n' % (each.strip()))
    handle.write('idba_ud --pre_correction --num_threads 1 --mink 20 --maxk 124 --step 20 --read %s --out %s_%s_k20-124\n' % (
    input_fasta, each_split[0], each_split[1]))
    handle.close()
