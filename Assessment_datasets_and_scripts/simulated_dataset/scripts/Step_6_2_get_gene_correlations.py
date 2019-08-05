import os
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-predicted_HGTs',
                    required=True,
                    help='MetaCHIP predicted HGTs (with direction)')

parser.add_argument('-t',
                    required=True,
                    help='Sequences of transferred genes in multi-fasta format')

parser.add_argument('-d',
                    required=True,
                    help='distribution of transfers to the recipient genomes')

args = vars(parser.parse_args())

HGT_candidates_seq = args['predicted_HGTs']
transfers_fasta = args['t']
distribution_of_transfers = args['d']

pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
pwd_blast_result = 'blast_result_get_gene_correlations.tab'
gene_correlations = 'HGT_candidates_transfer_correlations.txt'

# get transfer_to_recipient_dict
transfer_to_recipient_dict = {}
for each_recipient in open(distribution_of_transfers):
    each_recipient_split = each_recipient.strip().split(',')
    recipient_genome = each_recipient_split[0]
    transferred_genes = each_recipient_split[1:]
    for each_transferred_gene in transferred_genes:
        transfer_to_recipient_dict[each_transferred_gene] = recipient_genome

# run BlastN between MetaCHIP predicted candidates and introduced transfers
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, transfers_fasta))
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, HGT_candidates_seq, transfers_fasta, pwd_blast_result, outfmt))

# parse blast results
recovered_transfers = []
gene_correlations_handle = open(gene_correlations, 'w')
gene_correlations_handle.write('Recipient\tPredicted\tIntroduced\tIdentity\tAlignment_len\tPredicted_len\tIntroduced_len\n')
for match in open(pwd_blast_result):
    match_split = match.strip().split('\t')
    query = match_split[0]
    query_genome = query.split('_')[0]
    subject = match_split[1]
    subject_genome = subject.split('_')[0]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))
    if (identity >= 99) and (coverage_q >= 98) and (coverage_s >= 98):
        if (subject not in recovered_transfers) and ((query_genome == subject_genome) or (query_genome == transfer_to_recipient_dict[subject])):
            gene_correlations_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (transfer_to_recipient_dict[subject], query, subject, identity, align_len, query_len, subject_len))
            recovered_transfers.append(subject)
gene_correlations_handle.close()

# remove tmp files
os.system('rm *.nhr')
os.system('rm *.nin')
os.system('rm *.nog')
os.system('rm *.nsd')
os.system('rm *.nsi')
os.system('rm *.nsq')
print('The number of MetaCHIP recovered gene transfers: %i' % len(recovered_transfers))
