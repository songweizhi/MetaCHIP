import os
import argparse


def get_intersection(list_1, list_2):
    intersect_list = []
    for each in list_1:
        if (each in list_2) and (each not in intersect_list):
            intersect_list.append(each)
    return intersect_list


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

parser.add_argument('-mutation_level',
                    required=False,
                    type=int,
                    help='mutation level')

args = vars(parser.parse_args())
scaffold_file = args['a']
transfers_fasta = args['t']
minf = args['minf']
mutation_level = 'm%s' % args['mutation_level']


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
            if query not in recovered_transfers:
                recovered_transfers.append(query)
    if send < sstart:
        if (identity >= 99) and (coverage_q >= 99) and ((send >= minf) or (subject_len - sstart >= minf)):
            if query not in recovered_transfers:
                recovered_transfers.append(query)

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

############################## get all binned transfers ##############################

combined_ffn = 'combined_%s.ffn' % mutation_level
transfers_seq = transfers_fasta
blast_output = 'blast_results_binned.tab'
transfer_distribution = 'distribution_of_transfers.txt'

# run blast
os.system('makeblastdb -in %s -dbtype nucl -parse_seqids' % transfers_seq)
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
os.system('blastn -query %s -db %s -out %s %s' % (combined_ffn, transfers_seq, blast_output, blast_parameters))

transfer2recipient_dict = {}
for each_donor in open(transfer_distribution):
    each_donor_split = each_donor.strip().split(',')
    recipient_genome = each_donor_split[0]
    transfers = each_donor_split[1:]
    for each_transfer in transfers:
        transfer2recipient_dict[each_transfer] = recipient_genome

binned_transfers_overall = []
for each_hit in open(blast_output):
    each_hit_split = each_hit.strip().split('\t')
    query = each_hit_split[0]
    query_genome = query.split('_')[0]
    subject = each_hit_split[1]
    identity = float(each_hit_split[2])
    align_len = int(each_hit_split[3])
    query_len = int(each_hit_split[12])
    subject_len = int(each_hit_split[13])
    query_coverage = align_len/query_len
    subject_coverage = align_len/subject_len

    if (identity >= 99) and (query_coverage >= 0.99) and (subject_coverage > 0.99) and (transfer2recipient_dict[subject] == query_genome):
        if subject not in binned_transfers_overall:
            binned_transfers_overall.append(subject)


binned_transfers_from_recipient = get_intersection(binned_transfers_overall, recovered_transfers_m0)
binned_transfers_from_donor = get_intersection(binned_transfers_overall, recovered_transfers)
binned_transfers_from_both = get_intersection(binned_transfers_overall, recovered_transfers_both)

# report
print('Recovered transfers_donor: %s' % len(recovered_transfers_m0))
print('Recovered transfers_recipient: %s' % len(recovered_transfers))
print('Recovered transfers_intersection: %s' % len(recovered_transfers_both))
#print('%s,%s,%s' % (len(recovered_transfers_m0), len(recovered_transfers), len(recovered_transfers_both)))
#print('Recovered transfers(identity >= 99, coverage >= 99 and at least one flanking side longer than %sbp): %s' % (minf, len(recovered_transfers)))

print('binned_transfers_from_donor: %s' % len(binned_transfers_from_donor))
print('binned_transfers_from_recipient: %s' % len(binned_transfers_from_recipient))
print('binned_transfers_union: %s' % len(binned_transfers_overall))
print('binned_transfers_intersection: %s' % len(binned_transfers_from_both))

print('%s\t%s\t%s\t%s\t%s\t%s' % (len(recovered_transfers_m0),
                                  len(recovered_transfers),
                                  len(recovered_transfers_both),
                                  len(binned_transfers_from_donor),
                                  len(binned_transfers_from_recipient),
                                  len(binned_transfers_from_both)))
