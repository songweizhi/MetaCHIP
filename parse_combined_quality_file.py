import os

os.chdir('/Users/songweizhi/Desktop')

input_file = 'combined.txt'
bin_prefix = 'bin'

bin_id_2_lineage_dict = {}
for each in open(input_file):
    if each.startswith('  %s' % bin_prefix):
        #print(each.strip())
        each_split = each.strip().split(' ')
        #print(each_split)
        each_split_no_space = []
        for each_element in each_split:
            if each_element != '':
                each_split_no_space.append(each_element)
        #print(each_split_no_space)

        bin_id = each_split_no_space[0]
        lineage = each_split_no_space[1].split('__')[1]
        #print(lineage)
        bin_id_2_lineage_dict[bin_id] = lineage




for each_bin in open('grouping.txt'):
    each_bin_split = each_bin.strip().split(',')
    print('%s\t%s\t%s' % (each_bin_split[0], each_bin_split[1], bin_id_2_lineage_dict[each_bin_split[1]]))
