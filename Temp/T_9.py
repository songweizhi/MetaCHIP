
from Bio import SeqIO

# get id list of the 100 transfers

transfer_id_list = []
for each_transfer in SeqIO.parse('/Users/songweizhi/Desktop/sequences_of_gene_transfers.fasta', 'fasta'):
    transfer_id_list.append(each_transfer.id)


recovered_num = 0
for each in open('/Users/songweizhi/Desktop/HGT_candidates_with_direction.txt'):
    each_split = each.strip().split('\t')
    donor = each_split[0]
    recipient = each_split[1]
    if (donor in transfer_id_list) or (recipient in transfer_id_list):
        print(each.strip())
        recovered_num += 1

print(recovered_num)












