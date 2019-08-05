import os
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyfit


def get_group_size_and_hgt_num_list(genome_size_file, grouping_file, HGTs_PG_validated):
    genome_size_dict = {}
    for genome in open(genome_size_file):
        genome_split = genome.strip().split('\t')
        genome_name = '.'.join(genome_split[0].split('.')[:-1])
        if not genome.startswith('Genome'):
            genome_size_dict[genome_name] = float(genome_split[1])
    print(genome_size_dict)

    group_genome_dict = {}
    group_id_list = []
    for group in open(grouping_file):
        group_split = group.strip().split(',')
        group_id = group_split[0]
        genome_id = group_split[1]
        if group_id not in group_genome_dict:
            group_genome_dict[group_id] = [genome_id]
        else:
            group_genome_dict[group_id].append(genome_id)

        if group_id not in group_id_list:
            group_id_list.append(group_id)
    print(group_genome_dict)

    group_size_dict = {}
    for each_group in group_genome_dict:
        genome_list = group_genome_dict[each_group]
        total_len = 0
        for each_genome in genome_list:
            total_len += genome_size_dict[each_genome]
        group_size_dict[each_group] = total_len
    print(group_size_dict)

    group_hgt_num_dict = {}
    for hgt in open(HGTs_PG_validated):
        hgt_split = hgt.strip().split('\t')

        if hgt_split[2] not in group_hgt_num_dict:
            group_hgt_num_dict[hgt_split[2]] = 1
        else:
            group_hgt_num_dict[hgt_split[2]] += 1

        if hgt_split[3] not in group_hgt_num_dict:
            group_hgt_num_dict[hgt_split[3]] = 1
        else:
            group_hgt_num_dict[hgt_split[3]] += 1
    print(group_hgt_num_dict)

    group_id_list_sorted = sorted(group_id_list)

    group_hgt_num_list = []
    group_size_list = []
    for g in group_id_list_sorted:
        current_group_hgt_num = 0
        if g in group_hgt_num_dict:
            current_group_hgt_num = group_hgt_num_dict[g]

        group_hgt_num_list.append(current_group_hgt_num)
        group_size_list.append(group_size_dict[g])

    return group_id_list_sorted, group_size_list, group_hgt_num_list


wd = '/Users/songweizhi/Desktop/s16_sequence_len_HGT_num_correlation'
os.chdir(wd)


genome_size_file_human_gut = 'MetaBAT138bins_all_genome_size.txt'
genome_size_file_northsea = 'NorthSea37bins_all_genome_size.txt'

grouping_file_human_gut = 'MetaBAT138bins_o29_grouping.txt'
grouping_file_northsea = 'NorthSea37bins_o16_grouping.txt'

HGTs_PG_validated_human_gut = 'MetaBAT138bins_o29_HGTs_PG_validated.txt'
HGTs_PG_validated_northsea = 'NorthSea37bins_o16_HGTs_PG_validated.txt'


group_id_list_sorted_human_gut, group_size_list_human_gut, group_hgt_num_list_human_gut = get_group_size_and_hgt_num_list(genome_size_file_human_gut, grouping_file_human_gut, HGTs_PG_validated_human_gut)
group_id_list_sorted_northsea, group_size_list_northsea, group_hgt_num_list_northsea = get_group_size_and_hgt_num_list(genome_size_file_northsea, grouping_file_northsea, HGTs_PG_validated_northsea)


hgt_preference_human_gut = [num/size for num,size in zip(group_hgt_num_list_human_gut, group_size_list_human_gut)]
hgt_preference_northsea =  [num/size for num,size in zip(group_hgt_num_list_northsea, group_size_list_northsea)]


print(group_id_list_sorted_human_gut)
print(group_size_list_human_gut)
print(group_hgt_num_list_human_gut)
print(hgt_preference_human_gut)
print()
print(group_id_list_sorted_northsea)
print(group_size_list_northsea)
print(group_hgt_num_list_northsea)
print(hgt_preference_northsea)
print()


hgt_preference_human_gut_array = np.array(hgt_preference_human_gut)
hgt_preference_northsea_array = np.array(hgt_preference_northsea)

print(hgt_preference_human_gut_array)
print(hgt_preference_northsea_array)

k2_human_gut, p_human_gut = stats.normaltest(hgt_preference_human_gut_array)
k2_northsea, p_northsea = stats.normaltest(hgt_preference_northsea_array)

print(p_human_gut)  # not normal distribution
print(p_northsea)   # not normal distribution


print('####################')


x, y = stats.ranksums(hgt_preference_human_gut_array, hgt_preference_northsea_array)

print(x)
print(y)



print('####################')


x, y = stats.mannwhitneyu(hgt_preference_human_gut_array, hgt_preference_northsea_array)

print(x)
print(y)



# print(','.join(group_id_list_sorted))
# print(','.join(group_size_list_str))
# print(','.join(group_hgt_num_list_str))
#
#
# plt.scatter(group_size_list, group_hgt_num_list, label = group_id_list_sorted)
#
#
# for i, txt in enumerate(group_id_list_sorted):
#     plt.annotate(txt, (group_size_list[i], group_hgt_num_list[i]))
#
# plt.xlabel("Group size (Mbp)")
# plt.ylabel("HGT number")
#
# plt.show()
#



