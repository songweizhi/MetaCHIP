import os
import glob
import shutil

folders = ['AA', 'CF', 'ER', 'HO', 'PA', 'PC', 'RC', 'WT']

wd = '/srv/scratch/z5039045/3_fix_assemble_errors'
pwd_ra2 = '/srv/scratch/z5039045/Softwares/fix_assembly_errors/bin/ra2.py'
pwd_fix_fasta = '/srv/scratch/z5039045/Softwares/fix_assembly_errors/bin/fix_fasta.py'
pwd_nr_fasta = '/srv/scratch/z5039045/Softwares/fix_assembly_errors/bin/nr_fasta.py'
pwd_select_contig = '/srv/scratch/z5039045/Scripts/select_contig.pl'
pwd_get_fasta_stats = '/srv/scratch/z5039045/Scripts/get_fasta_stats.pl'
ra2_parameters = '-m 2 -c 2 --break -t 12'


pwd_qsub_file_folder = '%s/qsub_file_folder' % wd
if not os.path.exists(pwd_qsub_file_folder):
    os.makedirs(pwd_qsub_file_folder)
else:
    shutil.rmtree(pwd_qsub_file_folder)
    os.makedirs(pwd_qsub_file_folder)


nodes_number = 1
ppn_number = 12
memory = 120
walltime_needed = '23:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['python/3.4.3', 'shrinksam/0.9.0', 'velvet/1.2.10', 'bowtie/2.2.9']
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae\n\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'


for each in folders:
    pwd_each_folder = '%s/Refined_cf_bins_%s' % (wd, each)
    #print(pwd_each_folder)
    bins = '%s/*.fasta' % pwd_each_folder
    pwd_reads_file_R1 = '%s/Reads_files/%s_R1.fastq.gz' % (wd, each)
    pwd_reads_file_R2 = '%s/Reads_files/%s_R2.fastq.gz' % (wd, each)
    bin_names = [os.path.basename(file_name) for file_name in glob.glob(bins)]
    for each_bin in bin_names:
        #print(each_bin)
        each_bin_name, each_bin_ext = os.path.splitext(each_bin)

        pwd_qsub_file = '%s/%s.sh' % (pwd_qsub_file_folder, each_bin_name)
        qsub_file_handle = open(pwd_qsub_file, 'w')
        qsub_file_handle.write(header)
        qsub_file_handle.write(module_lines)
        #print(header)
        #print(module_lines)
        #print('cd %s\n' % wd)
        qsub_file_handle.write('\ncd %s\n' % wd)
        cmd_1 = 'python3 %s -i %s/%s -1 %s -2 %s %s\n' % (pwd_ra2, pwd_each_folder, each_bin, pwd_reads_file_R1, pwd_reads_file_R2, ra2_parameters)
        #print(cmd_1)
        qsub_file_handle.write(cmd_1)

        #print('cd %s/%s.curated' % (wd, each_bin_name))
        qsub_file_handle.write('cd %s/%s.curated\n' % (wd, each_bin_name))
        pwd_curated_bin_file = '%s/%s.curated/%s_curated%s' % (wd, each_bin_name, each_bin_name, each_bin_ext)
        pwd_curated_bin_file_lt2500 = '%s/%s.curated/%s_curated_lt2500%s' % (wd, each_bin_name, each_bin_name, each_bin_ext)
        pwd_curated_bin_file_lt2500_txt = '%s/%s.curated/%s_curated_lt2500.txt' % (wd, each_bin_name, each_bin_name)

        cmd_2 = 'python3 %s re_assembled.fa | python3 %s rename - > %s\n' % (pwd_fix_fasta, pwd_nr_fasta, pwd_curated_bin_file)
        #print(cmd_2)
        qsub_file_handle.write(cmd_2)
        cmd_3 = 'perl %s -m 2500 %s %s\n' % (pwd_select_contig, pwd_curated_bin_file, pwd_curated_bin_file_lt2500)
        #print(cmd_3)
        qsub_file_handle.write(cmd_3)
        cmd_4 = 'perl %s -T %s > %s\n' % (pwd_get_fasta_stats, pwd_curated_bin_file_lt2500, pwd_curated_bin_file_lt2500_txt)
        #print(cmd_4)
        qsub_file_handle.write(cmd_4)
        cmd_5 = 'perl %s -T %s/%s > %s/%s.curated/%s.txt\n' % (pwd_get_fasta_stats, pwd_each_folder, each_bin, wd, each_bin_name, each_bin_name)
        #print(cmd_5)
        qsub_file_handle.write(cmd_5)


        qsub_file_handle.close()

        #print('\n')
        os.chdir(pwd_qsub_file_folder)
        # os.system('qsub %s' % pwd_qsub_file)
        print('qsub %s' % pwd_qsub_file)



