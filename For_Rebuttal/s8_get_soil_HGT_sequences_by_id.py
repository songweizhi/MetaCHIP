import os
from Bio import SeqIO


hgt_fst_file = '/Users/songweizhi/Desktop/nature_rebuttal/nature10571_hgt_seqs/hgt.fst'
genome_meadata_file = '/Users/songweizhi/Desktop/nature_rebuttal/nature10571_genome_metadata.txt'
seqs_from_soil = '/Users/songweizhi/Desktop/nature_rebuttal/Soil_HGT_sequences_by_id.fasta'


seq_to_genome_dict = {}
for seq in SeqIO.parse(hgt_fst_file, 'fasta'):
    seq_description_split = str(seq.description).strip().split(' ')
    seq_to_genome_dict[seq_description_split[0]] = seq_description_split[1]


# get genome list
Soil_all_genome_list = set()
for genome_metadata in open(genome_meadata_file):
    genome_metadata_split = genome_metadata.strip().split('\t')
    genome_name = genome_metadata_split[0]
    niche = genome_metadata_split[2]
    if 'Soil' in niche:
        Soil_all_genome_list.add(genome_name)


Soil_all_sequence_list = set()
for seq_id in seq_to_genome_dict:
    if seq_to_genome_dict[seq_id] in Soil_all_genome_list:
        Soil_all_sequence_list.add(seq_id)


seq_num_all = 0
seq_num_soil = 0
seqs_from_soil_handle = open(seqs_from_soil, 'w')
for seq in SeqIO.parse(hgt_fst_file, 'fasta'):
    seq_id = str(seq.id)
    seq_num_all += 1
    if seq_id in Soil_all_sequence_list:
        seqs_from_soil_handle.write('>%s\n' % seq_id)
        seqs_from_soil_handle.write('%s\n' % str(seq.seq))
        seq_num_soil += 1
seqs_from_soil_handle.close()

print(seq_num_all)
print(seq_num_soil)
