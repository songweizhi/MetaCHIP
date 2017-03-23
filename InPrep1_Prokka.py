import os
import glob
import shutil
import argparse
import configparser

usage = """

    Usage:

    python3 path/to/this/script.py

    It will:
    1. Generate qsub file to run Prokka for each input genome
    2. Submit generated qsub files

"""

####################################################### Configuration ##################################################

parser = argparse.ArgumentParser()

parser.add_argument('-genome_folder',
                    required=True,
                    help='folder name of input genomes')

parser.add_argument('-genome_extension',
                    required=True,
                    help='genome_file_extension')

args = vars(parser.parse_args())
genome_folder = args['genome_folder']
genome_extension = args['genome_extension']

prokka_modules = ['perl/5.20.1', 'infernal/1.1.1', 'blast+/2.2.31', 'hmmer/3.1b2', 'prodigal/2.6.3', 'tbl2asn/25.3', 'prokka/1.12']

########################################################################################################################


# get working directory
working_directory = os.getcwd()

# get genome name list
genome_list = [os.path.basename(file_name) for file_name in glob.glob('%s/%s/*.%s' % (working_directory, genome_folder, genome_extension))]
print('%i genomes were found in total.' % len(genome_list))

if os.path.isdir('prokka_output'):
    shutil.rmtree('prokka_output')
    if os.path.isdir('prokka_output'):
        shutil.rmtree('prokka_output')
    os.mkdir('./prokka_output')
else:
    os.mkdir('./prokka_output')

# prepare prokka qsub file header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(1) + ':ppn=' + str(1) + '\n'
line_3 = '#PBS -l vmem=' + str(10) + 'gb\n'
line_4 = '#PBS -l walltime=' +'00:59:00' + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + 'weizhi.song@student.unsw.edu.au' + '\n'
line_7 = '#PBS -m ae'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

if os.path.isdir('qsub_files'):
    pass
else:
    os.mkdir('./qsub_files')

for genome in genome_list:
    genome_name = genome.split('.')[0]
    prokka_qsub_out = open('%s/qsub_files/qsub_prokka_%s.sh' % (working_directory, genome_name), 'w')
    # check whether folder exist
    if os.path.isdir('./prokka_output/%s' % genome_name):
        shutil.rmtree('./prokka_output/%s' % genome_name)
        os.mkdir('./prokka_output/%s' % genome_name)
    else:
        os.mkdir('./prokka_output/%s' % genome_name)

    prokka_qsub_out.write(header + '\n\n')
    for module in prokka_modules:
        prokka_qsub_out.write('module load %s\n' % module)
    prokka_qsub_out.write('prokka --force --metagenome --prefix %s --locustag %s --strain %s --outdir %s/prokka_output/%s %s/%s/%s'
                          % (genome_name, genome_name, genome_name, working_directory, genome_name, working_directory, genome_folder, genome))
    prokka_qsub_out.close()

# submit qsub files
qsub_file_list = [os.path.basename(file_name) for file_name in glob.glob('./qsub_files/*.sh')]
current_wd = os.getcwd()
os.chdir('./qsub_files')
for qsub_file in qsub_file_list:
    os.system('qsub %s' % qsub_file)
os.chdir(current_wd)
