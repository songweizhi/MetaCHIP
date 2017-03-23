import os
from Bio import SeqIO

# specify input files
wd = '/Users/songweizhi/Desktop/get_gene_correlations'
MetaCHIP_output = 'HGT_candidates_GemSIM.txt'
transfers_fasta = 'output_sequence_nc_90.fasta'
combined_ffn = 'combined_GemSIM.ffn'


GemSIM_hgts = open('%s/%s' % (wd, MetaCHIP_output))
pwd_transfers_fasta = '%s/%s' % (wd, transfers_fasta)
combined_ffn_GemSIM = '%s/%s' % (wd, combined_ffn)

# get the sequences of MetaCHIP predicted candidates
GemSIM_hgt_candidate_list = []
for each in GemSIM_hgts:
    each_split = each.strip().split('\t')
    for each_2 in each_split:
        if each_2 not in GemSIM_hgt_candidate_list:
            GemSIM_hgt_candidate_list.append(each_2)

GemSIM_candidate_fasta = '%s/GemSIM_candidate.fasta' % wd
GemSIM_handle = open(GemSIM_candidate_fasta, 'w')
for seq_record in SeqIO.parse(combined_ffn_GemSIM, 'fasta'):
    if seq_record.id in GemSIM_hgt_candidate_list:
        SeqIO.write(seq_record, GemSIM_handle, 'fasta')
GemSIM_handle.close()


# run BlastN between MetaCHIP predicted candidates and introduced transfers
pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
pwd_blast_result = '%s/blast_result.tab' % wd
# make blast database
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_transfers_fasta))
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
print('\nRunning Blast, be patient...')
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, GemSIM_candidate_fasta, pwd_transfers_fasta, pwd_blast_result, outfmt))


blast_results = open(pwd_blast_result)
n = 0
transfers_with_match = []
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
    if (identity >= 99) and (coverage_q >= 90) and (coverage_s >= 90):
        if subject not in transfers_with_match:
            transfers_with_match.append(subject)
        n += 1
print('transfers_with_match(%s):\n%s' % (len(transfers_with_match), transfers_with_match))
