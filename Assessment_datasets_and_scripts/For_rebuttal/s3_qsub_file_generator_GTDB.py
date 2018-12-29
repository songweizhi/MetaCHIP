import os
import shutil


##################################################### CONFIGURATION ####################################################

# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 16
memory = 120
walltime_needed = '11:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['perl/5.20.1', 'prodigal/2.6.3', 'pplacer/1.1.alpha16', 'hmmer/3.1b2', 'fasttree/2.1.10', 'gcc/6.2.0', 'gsl/2.1', 'fastani/1.1', 'R/3.4.2']
qsub_file_folder = '/Users/songweizhi/Desktop/qsub_files_GTDB'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes_renamed_sep'


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
n = 1
while n <= 50:
    current_qsub_file = '%s/qsub_GTDB_subfolder_%s.sh' % (qsub_file_folder, n)
    handle = open(current_qsub_file, 'w')
    handle.write('\n' + header + '\n')
    handle.write('module load python/2.7.10\ncd ~\n. mypythonenv_gtdbtk/bin/activate\n')
    handle.write(module_lines)
    handle.write('\ncd %s\n' % wd_on_katana)
    GTDB_cmd = 'gtdbtk classify_wf --cpus %s --genome_dir subfolder_%s --out_dir subfolder_%s_GTDB --extension fasta --prefix subfolder_%s' % (ppn_number, n, n, n)
    handle.write(GTDB_cmd)
    handle.close()
    n += 1

