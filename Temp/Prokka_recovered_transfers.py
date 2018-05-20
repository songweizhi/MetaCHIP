import os
import argparse


parser = argparse.ArgumentParser()


parser.add_argument('-a',
                    required=True,
                    help='Assemblies in multi-fasta format')

parser.add_argument('-t',
                    required=True,
                    help='Sequences of transferred genes in multi-fasta format')

args = vars(parser.parse_args())

scaffold_file = args['a']
transfers_fasta = args['t']

#wd = '/Users/songweizhi/Desktop/new_core'
#scaffold_file = 'scaffold_k20-124_lt2500.fa'
#transfers_fasta = 'output_sequence_nc_70.fasta'


# make blast database
wd = os.getcwd()
pwd_transfers_fasta = '%s/%s' % (wd, transfers_fasta)
pwd_scaffold_file = '%s/%s' % (wd, scaffold_file)
pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
pwd_blast_result = '%s/blast_result.tab' % wd
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_scaffold_file))
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
print('\nRunning Blast, be patient...')
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, pwd_transfers_fasta, pwd_scaffold_file, pwd_blast_result, outfmt))


blast_results = open(pwd_blast_result)
recovered_transfers = []
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
    if (identity == 100) and (coverage_q >= 98):
        #print(match.strip())
        if query not in recovered_transfers:
            recovered_transfers.append(query)

print('Recovered transfers(identity = 100 and coverage >= 98): %s' % len(recovered_transfers))
