import random

def get_common_stop_sequence():
    positive_strand_stop_codon = ['TAA', 'TAG', 'TGA']
    negative_strand_stop_codon = ['TTA', 'CTA', 'TCA']
    nc_bases = ['A', 'T', 'G', 'C']
    random_pos = random.sample(positive_strand_stop_codon, 3)
    random_neg = random.sample(negative_strand_stop_codon, 3)
    random_nc = random.sample(nc_bases, 4)
    common_stop_sequence = '%s%s%s%s%s%s%s%s%s%s' % (random_pos[0],
                                                     random_nc[0],
                                                     random_pos[1],
                                                     random_nc[1],
                                                     random_pos[2],
                                                     random_neg[0],
                                                     random_nc[2],
                                                     random_neg[1],
                                                     random_nc[3],
                                                     random_neg[2])
    return common_stop_sequence



a = get_common_stop_sequence()
print(a)
