import os
import shutil


##################################################### CONFIGURATION ####################################################

wd = '/Users/songweizhi/Desktop/nature_rebuttal/'
ncbi_assembly_summary = '/Users/songweizhi/Desktop/nature_rebuttal/assembly_summary.txt'
pwd_genome_metadata = '/Users/songweizhi/Desktop/nature_rebuttal/nature10571_genome_metadata.txt'
taxon_id_to_fna_file = 'taxon_id_to_fna.txt'


# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '00:59:00'
email = '244289990@qq.com'
modules_needed = ['python/2.7.12', 'blast+/2.6.0']
qsub_file_folder = '/Users/songweizhi/Desktop/qsub_files_download_genome'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes'

os.chdir(wd)


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


genome_name_set = set()
genome_taxon_id_set = set()
genome_name_taxon_dict = {}
for genome_meta in open(pwd_genome_metadata):
    if not genome_meta.startswith('Genome_Name'):
        genome_meta_split = genome_meta.strip().split('\t')
        genome_name = genome_meta_split[0]
        genome_taxon_id = genome_meta_split[1]
        genome_name_set.add(genome_name)
        genome_taxon_id_set.add(genome_taxon_id)
        if genome_taxon_id not in genome_name_taxon_dict:
            genome_name_taxon_dict[genome_taxon_id] = genome_name
print('Taxon ID number: %s' % len(genome_name_taxon_dict))


counted_id = set()
for each_entry in open(ncbi_assembly_summary):
    if not each_entry.startswith('#'):
        each_entry_split = each_entry.strip().split('\t')
        each_entry_taxon_id = each_entry_split[5]
        organism_name = each_entry_split[7]
        each_entry_ftp = each_entry_split[19]
        fna_gz_file_name = '%s_genomic.fna.gz' % each_entry_ftp.split('/')[-1]
        pwd_fna_gz = '%s/%s' % (each_entry_ftp, fna_gz_file_name)

        if (each_entry_taxon_id in genome_taxon_id_set) and ('metagenome' not in each_entry) and (each_entry_taxon_id not in counted_id):
            counted_id.add(each_entry_taxon_id)
            qsub_file = '%s/qsub_download_genome_%s.sh' % (qsub_file_folder, each_entry_ftp.split('/')[-1])
            qsub_file_handle = open(qsub_file, 'w')
            qsub_file_handle.write('\n' + header + '\n')
            qsub_file_handle.write(module_lines)
            qsub_file_handle.write('\ncd %s\n' % (wd_on_katana))
            download_cmd = 'wget %s' % pwd_fna_gz
            qsub_file_handle.write(download_cmd)
            qsub_file_handle.close()

            with open(taxon_id_to_fna_file, 'a') as file_handle:
                file_handle.write('%s\t%s\n' % (each_entry_taxon_id, fna_gz_file_name))
