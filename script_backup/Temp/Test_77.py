import os
import glob
import math
import random
import shutil
from Bio import SeqIO


wd = '/Users/songweizhi/Desktop/mNC'
subseq_num = 100
completeness = 0.3

os.chdir(wd)

# create output folder
output_folder = 'output_%s' % completeness
if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
    os.mkdir(output_folder)
else:
    os.mkdir(output_folder)


files = '*.fna'
file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]
print(file_list)

for input_file in file_list:
    # output file name
    file_name_no_extension, file_extension = os.path.splitext(input_file)
    output_file = '%s%s' % (file_name_no_extension, file_extension)


    # get the length of subsequences
    genome_length = 0
    genome_sequences = ''
    for each in SeqIO.parse('%s/%s' % (wd, input_file), 'fasta'):
        genome_length = len(each.seq)
        genome_sequences = str(each.seq)
    subseq_length = math.floor(genome_length/subseq_num)


    # get the positions to break the genome
    break_point_list = []
    m = 0
    while m < subseq_num:
        if m == subseq_num - 1:
            break_point_list.append([subseq_length*(subseq_num - 1) + 1, genome_length])
        else:
            break_point_list.append([subseq_length * m + 1, subseq_length * (m + 1)])
        m += 1


    # write out subsequences
    random_subseq_list = sorted(random.sample(break_point_list, int(subseq_num*completeness)))
    n = 1
    output_file_handle = open('%s/%s' % (output_folder, output_file), 'w')
    for each_subseq in random_subseq_list:
        current_subseq = genome_sequences[(each_subseq[0] - 1):(each_subseq[1] - 1)]
        current_subseq_id = 'strain_ctg%s' % n
        output_file_handle.write('>%s\n%s\n' % (current_subseq_id, current_subseq))
        n += 1
    output_file_handle.close()

