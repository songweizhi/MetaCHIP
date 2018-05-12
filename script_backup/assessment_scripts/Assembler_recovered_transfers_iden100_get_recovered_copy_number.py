import os

wd = '/Users/songweizhi/Desktop/assemblies'
introduced_transfers = 'input_sequence_mutant_nc_m0.fasta'
assemblies = 'm30_IDBA_UD_9_million_k20-124.fasta'
minf = 1000

pwd_blastn_exe = 'blastn'
pwd_makeblastdb_exe = 'makeblastdb'
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
pwd_blast_result = 'blast_result.tab'

os.chdir(wd)

os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, assemblies))
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, introduced_transfers, assemblies, pwd_blast_result, outfmt))

match_1 = []
match_2 = []
for match in open(pwd_blast_result):
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

    if (send < sstart) and (identity >= 99) and (coverage_q >= 99) and ((send >= minf) or (subject_len - sstart >= minf)):
        if query not in match_1:
            match_1.append(query)
        else:
            match_2.append(query)

    if (send > sstart) and (identity >= 99) and (coverage_q >= 99) and ((sstart >= minf) or (subject_len - send >= minf)):
        if query not in match_1:
            match_1.append(query)
        else:
            match_2.append(query)

print('\n')
print('Transfers recovered once: %s ' % len(match_1))
print(match_1)
print('Transfers recovered twice: %s' % len(match_2))
print(match_2)


