import os
__author__ = 'weizhisong'


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 12
memory = 80
walltime_needed = '3:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['bowtie/2.2.4']
output_folder = '/Users/weizhisong/Desktop/qsub_files/'

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
line_8 = '\n\ncd $PBS_O_WORDIR\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

################################################################

matches = open('/Users/weizhisong/Desktop/unitig_list.txt', 'r')
for unitig in matches:
    unitig = unitig.strip()
    unitig_name = unitig[:-4]
    part_1 = 'bowtie2 -x /srv/scratch/z5039045/PacBio/index_files/' + unitig_name
    part_2 = ' -f /srv/scratch/z3452659/ThomasPacBio-Nov15/data/2015-12-03.PacBioAssemblies/SON1053.SP16554.defaultsettings/SON1053.SP16554.hcq.preassembly.fasta'
    part_3 = ' -S /srv/scratch/z5039045/PacBio/sam_files/' + unitig_name
    output_file = open(output_folder + 'qsub_' + unitig_name + '.sh', 'w')
    output_file.write(header + module_lines + part_1 + part_2 + part_3)
    print('qdel ' + 'qsub_' + unitig_name + '.sh')
    output_file.close()
