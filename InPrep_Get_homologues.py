import os
import shutil
import argparse
import configparser

usage = """

    Usage:

    python3 path/to/this/script.py

    It will:
    1. Combine all gbk files (from Prokka) in defined folder
    2. Combine all ffn files (from Prokka) in defined folder
    3. Make blastdb with combined ffn file
    4. Generate qsub file to run all versus all blast of combined ffn
    5. Submit generated qsub file

"""

####################################################### Configuration ##################################################

path_to_makeblastdb_executable = 'makeblastdb'
blast_module = 'blast+/2.6.0'
prokka_output = 'prokka_output'
########################################################################################################################


# get combined gbk file
os.system('cat ./%s/*/*.gbk > combined.gbk' % prokka_output)

# get combined ffn file
os.system('cat ./%s/*/*.ffn > combined.ffn' % prokka_output)

# create blastdb folder and run makeblastdb
if not os.path.isdir('qsub_files'):
    os.mkdir('./qsub_files')

if os.path.isdir('blastdb'):
    shutil.rmtree('blastdb')
os.mkdir('blastdb')

os.system('cp combined.ffn ./blastdb/')
makeblastdb_cmd = path_to_makeblastdb_executable + ' -in ./blastdb/combined.ffn -dbtype nucl -parse_seqids'
os.system(makeblastdb_cmd)

# prepare blastn qsub file header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=1:ppn=3\n'
line_3 = '#PBS -l vmem=20gb\n'
line_4 = '#PBS -l walltime=11:59:00' + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M wythe1987@163.com\n'
line_7 = '#PBS -m ae\n'
line_8 = 'cd ' + os.getcwd()
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8

qsub_blastn_file_name = 'qsub_blastn_combined_ffn.sh'
blast_qsub_out = open('./qsub_files/' + qsub_blastn_file_name, 'w')
blast_qsub_out.write(header + '\n\n')
blast_qsub_out.write('module load %s\n' % blast_module)
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
blast_qsub_out.write('blastn -query %s/combined.ffn -db %s/blastdb/combined.ffn -out %s/all_vs_all_ffn.tab %s' % (os.getcwd(), os.getcwd(), os.getcwd(), blast_parameters))
blast_qsub_out.close()

#submit qsub files
os.system('qsub ./qsub_files/%s' % qsub_blastn_file_name)
