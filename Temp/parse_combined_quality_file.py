import os

os.chdir('/Users/songweizhi/Desktop')

input_file = 'combined.txt'
bin_prefix = 'bin'

bin_id_2_lineage_dict = {}
for each in open(input_file):
    if each.startswith('  %s' % bin_prefix):
        each_split = each.strip().split(' ')
        each_split_no_space = []
        for each_element in each_split:
            if each_element != '':
                each_split_no_space.append(each_element)

        bin_id = each_split_no_space[0]
        completeness = each_split_no_space[12]
        contamination = each_split_no_space[13]

        print('%s\t%s\t%s' % (bin_id, completeness, contamination))
