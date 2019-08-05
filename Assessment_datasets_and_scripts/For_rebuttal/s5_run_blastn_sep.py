import os
import shutil


##################################################### CONFIGURATION ####################################################

# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 3
memory = 20
walltime_needed = '11:59:00'
email = '244289990@qq.com'
modules_needed = ['blast+/2.6.0']
qsub_file_folder = '/Users/songweizhi/Desktop/qsub_files_blastn'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/Genomes2094_MetaCHIP_wd'


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
for taxon_id in open('/Users/songweizhi/Desktop/nature_rebuttal/taxon_id.txt'):

    taxon_id = taxon_id.strip()
    current_qsub_file = '%s/qsub_blastn_%s.sh' % (qsub_file_folder, taxon_id)
    handle = open(current_qsub_file, 'w')
    handle.write('\n' + header + '\n')
    handle.write(module_lines)
    handle.write('\ncd %s\n' % wd_on_katana)
    db_file = 'blastn_sep_wd/db_files/combined.ffn'
    blast_parameters = '-evalue 1e-5 -num_threads %s -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn' % ppn_number
    blastn_cmd = 'blastn -query Genomes2094_ffn_files/%s.ffn -db %s -out blastn_sep_wd/%s_blastn.tab %s' % (taxon_id, db_file, taxon_id, blast_parameters)
    handle.write(blastn_cmd)
    handle.close()


