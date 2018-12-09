import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['blast+/2.6.0']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_blastn_sep'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/blastn_sep'

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

################################################################
sample_prefix_file = 'ffn_files.txt'

for sample in open(sample_prefix_file):
    sample = sample.strip()
    blast_parameters = '-evalue 1e-5 -num_threads 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
    db = '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/GoodBins_0.5_0.05_blastdb/GoodBins_0.5_0.05_combined.ffn'
    cmd = 'blastn -query /srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/GoodBins_0.5_0.05_ffn_files/%s -db %s -out %s.blastn.tab %s' % (sample, db, sample, blast_parameters)

    print(cmd)



    output_handle = open('%s/qsub_blastn_sep_%s.sh' % (outputs_folder, sample), 'w')
    output_handle.write('%s\n' % header)
    output_handle.write(module_lines)
    output_handle.write('\ncd %s\n' % wd_on_katana)
    output_handle.write(cmd)
    output_handle.close()

