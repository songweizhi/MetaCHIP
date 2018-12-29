import re
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-predicted_HGTs',
                    required=True,
                    help='MetaCHIP predicted HGTs (with direction)')

parser.add_argument('-d',
                    required=True,
                    help='Distribution of transfers to the recipient genomes')

args = vars(parser.parse_args())
MetaCHIP_output = args['predicted_HGTs']
distribution_of_transfers = args['d']


# get transfer_to_recipient_dict
transferred_gene_id_list = []
transfer_to_recipient_dict = {}
for each_recipient in open(distribution_of_transfers):
    each_recipient_split = each_recipient.strip().split(',')
    recipient_genome = each_recipient_split[0]
    transferred_genes = each_recipient_split[1:]
    for each_transferred_gene in transferred_genes:
        transfer_to_recipient_dict[each_transferred_gene] = recipient_genome
        transferred_gene_id_list.append(each_transferred_gene)


detected_HGT = 0
for each_prediction in open(MetaCHIP_output):
    each_prediction_split = each_prediction.strip().split('\t')
    recipient = each_prediction_split[0]
    recipient_genome = '_'.join(recipient.split('_')[:-1])
    donor = each_prediction_split[1]
    donor_genome = '_'.join(donor.split('_')[:-1])
    # if (recipient in transferred_gene_id_list) or (donor in transferred_gene_id_list):
    #     detected_HGT += 1

    if recipient in transferred_gene_id_list:
        if transfer_to_recipient_dict[recipient] == donor_genome:
            detected_HGT += 1
    elif donor in transferred_gene_id_list:
        if transfer_to_recipient_dict[donor] == recipient_genome:
            detected_HGT += 1

print(detected_HGT)
