import os
from Bio import SeqIO

os.chdir('/Users/songweizhi/Desktop/tree_folder')
each_candidates = ['bin487_00383', 'bin292_01536']
pwd_blastp_exe = 'blastp'

HGT_genome_1 = '_'.join(each_candidates[0].split('_')[:-1])
HGT_genome_2 = '_'.join(each_candidates[1].split('_')[:-1])
blast_output = '%s___%s_gene_tree_blast.tab' % (each_candidates[0], each_candidates[1])
blast_output_sorted = '%s___%s_gene_tree_blast_sorted.tab' % (each_candidates[0], each_candidates[1])

# get the sequence of HGTs
gene_tree_seq = '%s___%s_gene_tree.seq' % (each_candidates[0], each_candidates[1])
gene_tree_seq_uniq = '%s___%s_gene_tree_uniq.seq' % (each_candidates[0], each_candidates[1])
self_seq = '%s___%s_gene_tree_selfseq.seq' % (each_candidates[0], each_candidates[1])
non_self_seq = '%s___%s_gene_tree_nonselfseq.seq' % (each_candidates[0], each_candidates[1])
self_seq_handle = open(self_seq, 'w')
non_self_seq_handle = open(non_self_seq, 'w')
for each_seq in SeqIO.parse(gene_tree_seq, 'fasta'):
    each_seq_genome_id = '_'.join(each_seq.id.split('_')[:-1])
    if each_seq.id in each_candidates:
        SeqIO.write(each_seq, self_seq_handle, 'fasta')
    elif (each_seq_genome_id != HGT_genome_1) and (each_seq_genome_id != HGT_genome_2):
        SeqIO.write(each_seq, non_self_seq_handle, 'fasta')
self_seq_handle.close()
non_self_seq_handle.close()

# run blast
os.system('%s -query %s -subject %s -outfmt 6 -out %s' % (pwd_blastp_exe, self_seq, non_self_seq, blast_output))
os.system('cat %s | sort > %s' % (blast_output, blast_output_sorted))

# get best match from each genome
current_query_subject_genome = ''
current_bit_score = 0
current_best_match = ''
best_match_list = []
for each_hit in open(blast_output_sorted):
    each_hit_split = each_hit.strip().split('\t')
    query = each_hit_split[0]
    subject = each_hit_split[1]
    subject_genome = '_'.join(subject.split('_')[:-1])
    query_subject_genome = '%s___%s' % (query, subject_genome)
    bit_score = float(each_hit_split[11])
    if current_query_subject_genome == '':
        current_query_subject_genome = query_subject_genome
        current_bit_score = bit_score
        current_best_match = subject
    elif current_query_subject_genome == query_subject_genome:
        if bit_score > current_bit_score:
            current_bit_score = bit_score
            current_best_match = subject
    elif current_query_subject_genome != query_subject_genome:
        best_match_list.append(current_best_match)
        current_query_subject_genome = query_subject_genome
        current_bit_score = bit_score
        current_best_match = subject
best_match_list.append(current_best_match)

# export sequences
gene_tree_seq_all = best_match_list + each_candidates
gene_tree_seq_uniq_handle = open(gene_tree_seq_uniq, 'w')
for each_seq2 in SeqIO.parse(gene_tree_seq, 'fasta'):
    if each_seq2.id in gene_tree_seq_all:
        SeqIO.write(each_seq2, gene_tree_seq_uniq_handle, 'fasta')
gene_tree_seq_uniq_handle.close()












