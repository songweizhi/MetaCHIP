import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-a',
                    required=True,
                    help='Assemblies in multi-fasta format')

parser.add_argument('-t',
                    required=True,
                    help='Sequences of transferred genes in multi-fasta format')

parser.add_argument('-minf',
                    required=False,
                    default=1000,
                    type=int,
                    help='minimum length of flanking sequences')

args = vars(parser.parse_args())
scaffold_file = args['a']
transfers_fasta = args['t']
minf = args['minf']

# make blast database
wd = os.getcwd()
pwd_transfers_fasta = '%s/%s' % (wd, transfers_fasta)
pwd_scaffold_file = '%s/%s' % (wd, scaffold_file)
pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'

transfers_fasta_m0 = 'input_sequence_mutant_nc_m0.fasta'
pwd_blast_result = '%s/blast_result.tab' % wd
pwd_blast_result_m0 = '%s/blast_result_m0.tab' % wd

os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_scaffold_file))
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
print('\nRunning Blast, be patient...')
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, pwd_transfers_fasta, pwd_scaffold_file, pwd_blast_result, outfmt))
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, transfers_fasta_m0, pwd_scaffold_file, pwd_blast_result_m0, outfmt))


blast_results = open(pwd_blast_result)
recovered_transfers = []
for match in blast_results:
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    sstart = int(match_split[8])
    send = int(match_split[9])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))
    if send > sstart:
        if (identity >= 99) and (coverage_q >= 99) and ((sstart >= minf) or (subject_len - send >= minf)):
            #print(match.strip())
            if query not in recovered_transfers:
                recovered_transfers.append(query)
    if send < sstart:
        if (identity >= 99) and (coverage_q >= 99) and ((send >= minf) or (subject_len - sstart >= minf)):
            #print(match.strip())
            if query not in recovered_transfers:
                recovered_transfers.append(query)

print(recovered_transfers)

recovered_transfers_m0 = []
recovered_transfers_both = []

for match_m0 in open(pwd_blast_result_m0):
    match_split_m0 = match_m0.strip().split('\t')
    query_m0 = match_split_m0[0]
    subject_m0 = match_split_m0[1]
    identity_m0 = float(match_split_m0[2])
    align_len_m0 = int(match_split_m0[3])
    sstart_m0 = int(match_split_m0[8])
    send_m0 = int(match_split_m0[9])
    query_len_m0 = int(match_split_m0[12])
    subject_len_m0 = int(match_split_m0[13])
    coverage_q_m0 = float("{0:.2f}".format(float(align_len_m0) * 100 / float(query_len_m0)))
    coverage_s_m0 = float("{0:.2f}".format(float(align_len_m0) * 100 / float(subject_len_m0)))
    if send_m0 > sstart_m0:
        if (identity_m0 >= 99) and (coverage_q_m0 >= 99) and ((sstart_m0 >= minf) or (subject_len_m0 - send_m0 >= minf)):
            if query_m0 not in recovered_transfers_m0:
                recovered_transfers_m0.append(query_m0)
            if (query_m0 in recovered_transfers) and (query_m0 not in recovered_transfers_both):
                recovered_transfers_both.append(query_m0)

    if send_m0 < sstart_m0:
        if (identity_m0 >= 99) and (coverage_q_m0 >= 99) and ((send_m0 >= minf) or (subject_len_m0 - sstart_m0 >= minf)):
            if query_m0 not in recovered_transfers_m0:
                recovered_transfers_m0.append(query_m0)
            if (query_m0 in recovered_transfers) and (query_m0 not in recovered_transfers_both):
                recovered_transfers_both.append(query_m0)




print('Recovered transfers_m0: %s' % len(recovered_transfers_m0))
print('Recovered transfers: %s' % len(recovered_transfers))
print('Recovered transfers_both: %s' % len(recovered_transfers_both))


print('%s,%s,%s' % (len(recovered_transfers_m0), len(recovered_transfers), len(recovered_transfers_both)))

#print('Recovered transfers(identity >= 99, coverage >= 99 and at least one flanking side longer than %sbp): %s' % (minf, len(recovered_transfers)))



