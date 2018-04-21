
out_handle = open('/Users/songweizhi/Desktop/grouping_1_no_index.txt', 'w')

for each in open('/Users/songweizhi/Desktop/grouping_1.txt'):
    each_split = each.strip().split(',')
    group_index = each_split[0]
    group_id = group_index.split('_')[0]
    bin_id = each_split[1]
    #print('%s,%s' % (group_id, bin_id))
    out_handle.write('%s,%s\n' % (group_id, bin_id))
out_handle.close()


grouping_file = '/Users/songweizhi/Desktop/grouping_1_no_index.txt'



def get_number_of_group(grouping_file):

    group_list = []
    for each_genome in open(grouping_file):
        each_genome_split = each_genome.strip().split(',')
        group_id = each_genome_split[0]
        if group_id not in group_list:
            group_list.append(group_id)
    number_of_group = len(group_list)

    return number_of_group



print(get_number_of_group(grouping_file))