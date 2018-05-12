import os
from Bio import SeqIO


GemSIM_hgts = open('/Users/songweizhi/Desktop/get_gene_correlations/HGT_candidates_GemSIM.txt')
PC_hgts = open('/Users/songweizhi/Desktop/get_gene_correlations/HGT_candidates_PC.txt')


GemSIM_hgt_candidate_list = []
for each in GemSIM_hgts:
    each_split = each.strip().split('\t')
    for each_2 in each_split:
        if each_2 not in GemSIM_hgt_candidate_list:
            GemSIM_hgt_candidate_list.append(each_2)
print(GemSIM_hgt_candidate_list)
print(len(GemSIM_hgt_candidate_list))

PC_hgt_candidate_list = []
for each in PC_hgts:
    each_split = each.strip().split('\t')
    for each_2 in each_split:
        if each_2 not in PC_hgt_candidate_list:
            PC_hgt_candidate_list.append(each_2)
print(PC_hgt_candidate_list)
print(len(PC_hgt_candidate_list))


PC_candidate_fasta = '/Users/songweizhi/Desktop/get_gene_correlations/PC_candidate.fasta'
PC_handle = open(PC_candidate_fasta, 'w')
for seq_record in SeqIO.parse('/Users/songweizhi/Desktop/get_gene_correlations/combined_PC.ffn', 'fasta'):
    if seq_record.id in PC_hgt_candidate_list:
        SeqIO.write(seq_record, PC_handle, 'fasta')
PC_handle.close()


GemSIM_candidate_fasta = '/Users/songweizhi/Desktop/get_gene_correlations/GemSIM_candidate.fasta'
GemSIM_handle = open(GemSIM_candidate_fasta, 'w')
for seq_record in SeqIO.parse('/Users/songweizhi/Desktop/get_gene_correlations/combined_GemSIM.ffn', 'fasta'):
    if seq_record.id in GemSIM_hgt_candidate_list:
        SeqIO.write(seq_record, GemSIM_handle, 'fasta')
GemSIM_handle.close()


# run blast
pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
pwd_blast_result = '/Users/songweizhi/Desktop/get_gene_correlations/blast_result.tab'
# make blast database
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, PC_candidate_fasta))
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
print('\nRunning Blast, be patient...')
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, GemSIM_candidate_fasta, PC_candidate_fasta, pwd_blast_result, outfmt))


blast_results = open(pwd_blast_result)
n = 0
for match in blast_results:
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))

    if (identity >= 98) and (coverage_q >= 90) and (coverage_s >= 90):
        print(match.strip())
        n += 1
print(n)














