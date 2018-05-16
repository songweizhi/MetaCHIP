import os


wd  = '/Users/songweizhi/Desktop/PC_group'
nc_transfers = 'HGT_candidates_ET_mNC_PC.txt'
introduced_distribution = 'distribution_of_transfers.txt'

mutation_groups = [0, 5, 10, 15, 20, 25, 30]

os.chdir(wd)


# get already existed gene transfers
original_transfer_list = []
for each_original in open(nc_transfers):
    recipient = each_original.strip().split('\t')[0]
    donor = each_original.strip().split('\t')[1]
    if sorted([recipient, donor])[0] not in original_transfer_list:
        original_transfer_list.append(sorted([recipient, donor])[0])


# get introduced gene list
introduced_hgt_list = []
for each_recipient_genome in open(introduced_distribution):
    each_recipient_genome_split = each_recipient_genome.strip().split(',')
    transfered_genes = each_recipient_genome_split[1:]
    for each_transfer in transfered_genes:
        introduced_hgt_list.append(each_transfer)


print('Group\tTotal\tOriginal\tPredicted_original\tIntroduced\tPredicted_introduced\tNew_intriduced')
for each_mutation_group in mutation_groups:
    file_name = 'HGT_candidates_ET_m%s_PC.txt' % each_mutation_group
    total = []
    predicted_original_transfers = []
    predicted_introduced_transfers = []
    new_intriduced_transfers = []
    for each_prediction in open(file_name):
        each_prediction_recipient = each_prediction.strip().split('\t')[0]
        each_prediction_donor = each_prediction.strip().split('\t')[1]
        A = sorted([each_prediction_recipient, each_prediction_donor])[0]
        if (A in original_transfer_list) and (A not in predicted_original_transfers):
            predicted_original_transfers.append(A)
        if (A in introduced_hgt_list) and (A not in predicted_introduced_transfers):
            predicted_introduced_transfers.append(A)
        if (A not in original_transfer_list) and (A not in introduced_hgt_list):
            new_intriduced_transfers.append(A)
        if A not in total:
            total.append(A)

    print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_mutation_group,
                                          len(total),
                                          len(original_transfer_list),
                                          len(predicted_original_transfers),
                                          len(introduced_hgt_list),
                                          len(predicted_introduced_transfers),
                                          len(new_intriduced_transfers)))















