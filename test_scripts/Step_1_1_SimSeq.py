import os
import shutil
from Bio import SeqIO


def get_genome_size(fasta_file):
    genome = SeqIO.parse(fasta_file, 'fasta')
    total_length = 0
    for each_contig in genome:
        total_length += len(each_contig.seq)
    return total_length


wd = '/Users/songweizhi/Desktop/SimSeq_wd'
genome_folder = 'input_genomes'
abundance_file = 'abundance.txt'
reads_number = 4000000
SimSeq_wd = '1_SimSeq'
pwd_SimSeq_jar = '/srv/scratch/z5039045/Softwares/SimSeq/jars/SimSeq.jar'
pwd_error_model = '/srv/scratch/z5039045/Softwares/SimSeq-master/profiles/miseq_250bp.txt'
reads_length = 250
insert_size = 200
nodes_number = 1
ppn_number = 1
memory = 30
walltime_needed = '00:11:59'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['samtools/1.3.1', 'java/8u91']

pwd_genome_folder = '%s/%s' % (wd, genome_folder)
pwd_abundance_file = '%s/%s' % (wd, abundance_file)
pwd_qsub_file_folder = '%s/qsub_files' % wd
pwd_SimSeq_wd = '%s/%s' % (wd, SimSeq_wd)

###################################### Prepare qsub file header ######################################

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae\n\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
#modules_needed_split = modules_needed.split(',')
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

######################################## Prepare input files ########################################


# get genome size list
abundance_combined = open(pwd_abundance_file)
genome_name_list = []
genome_name_no_extension_list = []
genome_size_list = []
replicate_number = 0
for each in open(pwd_abundance_file):
    each_split = each.strip().split('\t')
    genome_name = each_split[0]
    genome_name_no_extension = genome_name.split('.')[0]
    genome_file = '%s/%s/%s' % (wd, genome_folder, genome_name)
    genome_size = get_genome_size(genome_file)
    genome_size_mbp = float("{0:.2f}".format(genome_size / (1024 * 1024)))

    replicate_number = len(each_split) - 1
    genome_name_no_extension_list.append(genome_name_no_extension)
    genome_name_list.append(genome_name)

    genome_size_list.append(genome_size_mbp)

# prepare GemSIM working directory
if os.path.isdir(pwd_SimSeq_wd):
    shutil.rmtree(pwd_SimSeq_wd)
    os.mkdir(pwd_SimSeq_wd)
else:
    os.mkdir(pwd_SimSeq_wd)

if os.path.isdir(pwd_qsub_file_folder):
    shutil.rmtree(pwd_qsub_file_folder)
    os.mkdir(pwd_qsub_file_folder)
else:
    os.mkdir(pwd_qsub_file_folder)


n = 1
while n <= replicate_number:
    current_abundance_list = []
    for each in open(pwd_abundance_file):
        each_split = each.strip().split('\t')
        current_abundance = int(each_split[n])
        current_abundance_list.append(current_abundance)
    print(genome_name_no_extension_list)

    m = 0
    weight_list = []
    while m < len(genome_name_no_extension_list):
        current_weight = (genome_size_list[m])*(current_abundance_list[m])
        weight_list.append(current_weight)
        m += 1
    weight_summation = sum(weight_list)
    multiply_unit = reads_number/weight_summation
    reads_number_list = []
    for each in weight_list:
        current_reads_number = multiply_unit * each
        reads_number_list.append(int(current_reads_number))
    print(reads_number_list)

    pwd_current_replicate_wd = '%s/replicate%s_wd' % (pwd_SimSeq_wd, n)
    os.mkdir(pwd_current_replicate_wd)

    p = 0
    for current_genome in genome_name_list:
        current_genome_no_extension = current_genome.split('.')[0]
        pwd_current_genome = '%s/%s/%s' % (wd, genome_folder, current_genome)
        current_qsub_file = '%s/%s_replicate%s.sh' % (pwd_qsub_file_folder, current_genome_no_extension, n)
        current_qsub_file_handle = open(current_qsub_file, 'w')

        print('%s\t%s' % (current_genome, reads_number_list[p]))

        current_qsub_file_handle.write(header + module_lines)
        current_qsub_file_handle.write('cd %s\n' % pwd_current_replicate_wd)

        current_qsub_file_handle.write('java -jar -Xmx2048m %s -1 %s -2 %s -e %s -n %s -l %s -o %s.sam -r %s -p %s_\n' % (pwd_SimSeq_jar,
                                                                                               reads_length,
                                                                                               reads_length,
                                                                                               pwd_error_model,
                                                                                               reads_number_list[p],
                                                                                               insert_size,
                                                                                               current_genome_no_extension,
                                                                                               pwd_current_genome,
                                                                                               current_genome_no_extension))

        current_qsub_file_handle.write('samtools view -bS -T %s -o %s.bam %s.sam\n' % (pwd_current_genome,
                                                            current_genome_no_extension,
                                                            current_genome_no_extension))

        current_qsub_file_handle.write('samtools sort -o %s.sorted.bam %s.bam\n' % (current_genome_no_extension,
                                                         current_genome_no_extension))

        current_qsub_file_handle.write('samtools index %s.sorted.bam\n' % current_genome_no_extension)

        current_qsub_file_handle.write('samtools fastq -1 %s_R1.fastq -2 %s_R2.fastq %s.sorted.bam --reference %s\n' % (current_genome_no_extension,
                                                                                             current_genome_no_extension,
                                                                                             current_genome_no_extension,
                                                                                             pwd_current_genome))
        current_qsub_file_handle.close()
        current_wd = os.getcwd()
        os.chdir(pwd_qsub_file_folder)
        os.system('qsub %s' % current_qsub_file)
        os.chdir(current_wd)

        p += 1

    n += 1




