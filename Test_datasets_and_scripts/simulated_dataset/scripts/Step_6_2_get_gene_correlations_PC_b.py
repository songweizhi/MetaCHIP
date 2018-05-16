import re
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-predicted_HGTs',
                    required=True,
                    help='MetaCHIP predicted HGTs (with direction)')

parser.add_argument('-d',
                    required=True,
                    help='distribution of transfers to the recipient genomes')

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


blast_recovered = 0
blast_recovered_right_direction = 0
tree_validated = 0
tree_validated_right_direction = 0
both_in_right_direction = 0
for each_prediction in open(MetaCHIP_output):
    each_prediction_split = each_prediction.strip().split('\t')
    recipient = each_prediction_split[0]
    recipient_genome = '_'.join(recipient.split('_')[:-1])
    donor = each_prediction_split[1]
    donor_genome = '_'.join(donor.split('_')[:-1])
    direction_Blast = each_prediction_split[6]
    direction_Tree = each_prediction_split[7]
    if (recipient in transferred_gene_id_list) or (donor in transferred_gene_id_list):
        #print(each_prediction.strip())
        blast_recovered += 1
        if re.match('B[A-Z]*<-A[A-Z]*', direction_Blast):
            blast_recovered_right_direction += 1
        if direction_Tree != 'N/A':
            tree_validated += 1
            if re.match('A[A-Z]*-->B[A-Z]*', direction_Tree):
                tree_validated_right_direction += 1
        if (re.match('B[A-Z]*<-A[A-Z]*', direction_Blast)) and (re.match('A[A-Z]*-->B[A-Z]*', direction_Tree)):
            both_in_right_direction += 1

# print('blast_recovered: %s' % blast_recovered)
# print('blast_recovered_right_direction: %s' % blast_recovered_right_direction)
# print('tree_validated: %s' % tree_validated)
# print('tree_validated_right_direction: %s' % tree_validated_right_direction)

print(blast_recovered)
print(blast_recovered_right_direction)
print(tree_validated)
print(tree_validated_right_direction)
print(both_in_right_direction)
