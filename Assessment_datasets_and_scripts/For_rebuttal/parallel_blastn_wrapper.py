import os
import glob
import shutil
import argparse


usage = '''

module load python/3.5.2
python /srv/scratch/z5039045/Softwares/MetaCHIP/parallel_blastn_wrapper.py -p NorthSea -qsub 
cat Total2094_all_parallel_blastn_output/*.tab > Total2094_all_vs_all_ffn.tab

'''


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()
parser.add_argument('-p', required=True, help='prefix')
parser.add_argument('-qsub', required=False, action='store_true', help='prefix')

args = vars(parser.parse_args())
prefix = args['p']
submit_jobs = args['qsub']

nodes_number = 1
ppn_number = 3
memory = 20
walltime_needed = '11:59:00'
modules_needed = ['blast+/2.6.0']


################################################### define file name ###################################################

blastn_db = '%s_all_blastdb/%s_all_combined_ffn.fasta' % (prefix, prefix)
ffn_file_folder = '%s_all_prodigal_output' % (prefix)
qsub_file_folder = '%s_all_parallel_blastn_qsubs' % (prefix)
blast_results_folder = '%s_all_parallel_blastn_output' % (prefix)

# Create folder
force_create_folder(qsub_file_folder)
force_create_folder(blast_results_folder)

# get ffn file list
ffn_file_re = '%s/*.ffn' % ffn_file_folder
ffn_file_list = [os.path.basename(file_name) for file_name in glob.glob(ffn_file_re)]


################################################## generate_qsub_file ##################################################

# Prepare header
line_1 = '#!/bin/bash'
line_2 = '#PBS -l nodes=%s:ppn=%s' % (str(nodes_number), str(ppn_number))
line_3 = '#PBS -l vmem=%sgb' % str(memory)
line_4 = '#PBS -l walltime=%s' % walltime_needed
line_5 = '#PBS -j oe'
line_6 = '#PBS -m ae'
header = '%s\n%s\n%s\n%s\n%s\n%s\n' % (line_1, line_2, line_3, line_4, line_5, line_6)

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load %s\n' % module

# write to qsub files
for ffn_file in ffn_file_list:

    pwd_ffn_file = '%s/%s' % (ffn_file_folder, ffn_file)
    ffn_file_basename = '.'.join(ffn_file.split('.')[:-1])
    current_qsub_file = '%s/qsub_blastn_%s.sh' % (qsub_file_folder, ffn_file_basename)
    handle = open(current_qsub_file, 'w')
    handle.write('\n' + header + '\n')
    handle.write(module_lines)
    handle.write('\ncd %s\n' % os.getcwd())
    blast_parameters = '-evalue 1e-5 -num_threads %s -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn' % ppn_number
    blastn_cmd = 'blastn -query %s -db %s -out %s/%s_blastn.tab %s\n' % (pwd_ffn_file, blastn_db, blast_results_folder, ffn_file_basename, blast_parameters)
    handle.write(blastn_cmd)
    handle.close()

    # submit job
    if submit_jobs is True:
        os.system('qsub %s' % current_qsub_file)


######################################################## report ########################################################

cat_cmd = 'cat %s/*.tab > %s_all_vs_all_ffn.tab' % (blast_results_folder, prefix)
report_message = 'Combine blast results after all submitted jobs were finished with: \n%s' % cat_cmd

print(report_message)

