import os

os.chdir('/Users/songweizhi/Desktop/get_successfully_binned_transfers')
mutation_level = 'm30'
# codes for m0 need to be modified

combined_ffn = 'combined_%s.ffn' % mutation_level
transfers_seq = 'input_sequence_mutant_nc_%s.fasta' % mutation_level
blast_output = 'blast_results_binned_%s.tab' % mutation_level
transfer_distribution = 'distribution_of_transfers.txt'
assemblies = '%s_IDBA_UD_9_million_k20-124.fasta' % mutation_level
minf = 1000

############################## get all binned transfers ##############################

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

    if mutation_level == 'm0':
        if (identity >= 99) and (query_coverage >= 0.99) and (subject_coverage > 0.99) and (transfer2recipient_dict[subject] == query_genome):
            if subject not in binned_transfers_overall:
                binned_transfers_overall.append(subject)
            print(each_hit.strip())

    else:
        if (identity >= 99) and (query_coverage >= 0.99) and (subject_coverage > 0.99) and (transfer2recipient_dict[subject] == query_genome):
            if subject not in binned_transfers_overall:
                binned_transfers_overall.append(subject)

print(len(binned_transfers_overall))
print(binned_transfers_overall)




