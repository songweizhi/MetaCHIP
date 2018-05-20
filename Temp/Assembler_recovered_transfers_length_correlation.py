import os
import argparse
from Bio import SeqIO

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


########## get length distribution ##########
list_lt600 = []
list_lt800 = []
list_lt1000 = []
list_lt1200 = []
list_lt1400 = []
list_lt1600 = []
list_lt1800 = []
list_lt2000 = []
for each in SeqIO.parse(pwd_transfers_fasta, 'fasta'):
    seq_id = each.id
    seq_length = len(each.seq)
    if 600 <=seq_length < 800:
        list_lt600.append(seq_id)
    elif 800 <=seq_length < 1000:
        list_lt800.append(seq_id)
    elif 1000 <=seq_length < 1200:
        list_lt1000.append(seq_id)
    elif 1200 <=seq_length < 1400:
        list_lt1200.append(seq_id)
    elif 1400 <=seq_length < 1600:
        list_lt1400.append(seq_id)
    elif 1600 <=seq_length < 1800:
        list_lt1600.append(seq_id)
    elif 1800 <=seq_length < 2000:
        list_lt1800.append(seq_id)
    elif seq_length >= 2000:
        list_lt2000.append(seq_id)

# print('list_lt600: %s' % list_lt600)
# print('list_lt800: %s' % list_lt800)
# print('list_lt1000: %s' % list_lt1000)
# print('list_lt1200: %s' % list_lt1200)
# print('list_lt1400: %s' % list_lt1400)
# print('list_lt1600: %s' % list_lt1600)
# print('list_lt1800: %s' % list_lt1800)
# print('list_lt2000: %s' % list_lt2000)

# get recovered genes for each group
list_lt600_recovered = []
list_lt800_recovered = []
list_lt1000_recovered = []
list_lt1200_recovered = []
list_lt1400_recovered = []
list_lt1600_recovered = []
list_lt1800_recovered = []
list_lt2000_recovered = []
print(recovered_transfers)
print(len(recovered_transfers))
for each in recovered_transfers:
    if each in list_lt600:
        list_lt600_recovered.append(each)
    elif each in list_lt800:
        list_lt800_recovered.append(each)
    elif each in list_lt1000:
        list_lt1000_recovered.append(each)
    elif each in list_lt1200:
        list_lt1200_recovered.append(each)
    elif each in list_lt1400:
        list_lt1400_recovered.append(each)
    elif each in list_lt1600:
        list_lt1600_recovered.append(each)
    elif each in list_lt1800:
        list_lt1800_recovered.append(each)
    elif each in list_lt2000:
        list_lt2000_recovered.append(each)

print('600-800bp: %s' % (round(len(list_lt600_recovered)/len(list_lt600), 2)))
print('800-1000bp: %s' % (round(len(list_lt800_recovered)/len(list_lt800), 2)))
print('1000-1200bp: %s' % (round(len(list_lt1000_recovered)/len(list_lt1000), 2)))
print('1200-1400bp: %s' % (round(len(list_lt1200_recovered)/len(list_lt1200), 2)))
print('1400-1600bp: %s' % (round(len(list_lt1400_recovered)/len(list_lt1400), 2)))
print('1600-1800bp: %s' % (round(len(list_lt1600_recovered)/len(list_lt1600), 2)))
print('1800-2000bp: %s' % (round(len(list_lt1800_recovered)/len(list_lt1800), 2)))
print('>2000bp: %s' % (round(len(list_lt2000_recovered)/len(list_lt2000), 2)))

#print('Recovered transfers(identity >= 99, coverage >= 99 and at least one flanking side longer than %sbp): %s' % (minf, len(recovered_transfers)))









