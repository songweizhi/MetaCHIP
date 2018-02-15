import os
__author__ = 'weizhisong'


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 20
walltime_needed = '11:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['python/3.5.2', 'perl/5.20.1', 'blast+/2.6.0']
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

for each in open('/Users/songweizhi/Desktop/47bins.txt'):
    each = each.strip()
    handle = open('/Users/songweizhi/Desktop/qsub_files/qsub_COG_%s.sh' % each, 'w')
    handle.write('\n' + header + '\n')
    handle.write(module_lines)
    handle.write('cd /srv/scratch/z5039045/MetaCHIP/BBAY68_69_70/47bins_MetaCHIP_wd/prokka_output/%s\n' % each)
    handle.write('python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in %s.faa -t P\n' % each)
    handle.close()
