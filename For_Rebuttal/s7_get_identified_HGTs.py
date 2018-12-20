
import os
from Bio import SeqIO
import matplotlib.pyplot as plt


wd = '/Users/songweizhi/Desktop/check'
os.chdir(wd)


# input files
query_seq = 'Soil368Genomes_g139_HGTs_PG_nc.fasta'
HGT_candidates_BM = 'Soil368Genomes_g139_HGTs_PG.txt'
HGT_candidates_PG = 'Soil368Genomes_g139_HGTs_PG_validated.txt'


subject_seq = '/Users/songweizhi/Desktop/nature_rebuttal/Soil_HGT_sequences_by_id.fasta'

# define file name
recent_HGT_seq = 'recent_HGT.fasta'
recent_hgt_blast_output = 'output_recent_hgt.tab'
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'


# get MetaCHIP identified recent transfers
recent_HGT_num = 0
recent_hgts = set()
identity_list_all = []
for each in open(HGT_candidates_PG):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        identity = float(each_split[4])
        identity_list_all.append(identity)
        if identity >= 99:
            recent_HGT_num += 1
            recent_hgts.add(gene_1)
            recent_hgts.add(gene_2)
print('The number of MetaCHIP identified recent HGTs: %s' % recent_HGT_num)
print('The number of genes involving recent HGT: %s' % len(recent_hgts))


# get sequence of MetaCHIP identified recent transfers
recent_HGT_seq_handle = open(recent_HGT_seq, 'w')
for each in SeqIO.parse(query_seq, 'fasta'):
    if str(each.id) in recent_hgts:
        recent_HGT_seq_handle.write('>%s\n' % str(each.id))
        recent_HGT_seq_handle.write('%s\n' % str(each.seq))
recent_HGT_seq_handle.close()


# run blast between MetaCHIP identified recent transfers and pubished results
recent_hgt_blast_command = 'blastn -query %s -subject %s -out %s %s' % (recent_HGT_seq, subject_seq, recent_hgt_blast_output, blast_parameters)
os.system(recent_hgt_blast_command)


# parse blast results for recent HGTs (should blast genes from recipients only, current both donor and recipient genes)
validated_recent_hgt = set()
for blast_hit in open(recent_hgt_blast_output):
    blast_hit_split = blast_hit.strip().split('\t')
    query = blast_hit_split[0]
    identity = float(blast_hit_split[2])
    align_len = int(blast_hit_split[3])
    query_len = int(blast_hit_split[12])
    subject_len = int(blast_hit_split[13])
    query_cov = align_len / query_len
    subject_cov = align_len / subject_len
    if (identity == 100) and (query_cov == 1):
        validated_recent_hgt.add(query)


print(len(validated_recent_hgt))


# plot identity distribution
# plt.hist(identity_list_all, bins=30, linewidth=1, alpha=0.5)
# plt.show()



