
pwd_HGT_candidates_file = '/Users/songweizhi/Desktop/HGT_candidates.txt'

pwd_HGT_candidates_with_direction_file = '/Users/songweizhi/Desktop/HGT_candidates_with_direction.txt'

def add_direction(input_file, output_file):
    input = open(input_file)
    output = open(output_file, 'w')

    # get overall list
    overall = []
    for each in input:
        overall.append(each.strip())

    # get overlap list
    tmp_list = []
    overlap_list = []
    for each in overall:
        each_split = each.split('\t')
        each_reverse = '%s\t%s' % (each_split[1], each_split[0])
        tmp_list.append(each)
        if each_reverse in tmp_list:
            overlap_list.append(each)

    # get non-overlap list
    non_overlap_list = []
    for each in overall:
        each_split = each.split('\t')
        each_reverse = '%s\t%s' % (each_split[1], each_split[0])
        if (each not in overlap_list) and (each_reverse not in overlap_list):
            non_overlap_list.append(each)

    # get output
    for each in non_overlap_list:
        each_split = each.split('\t')
        each_split_0 = each.split('\t')[0]
        each_split_1 = each.split('\t')[1]
        each_split_0_bin = each_split_0.split('_')[0]
        each_split_1_bin = each_split_1.split('_')[0]
        output.write('%s\t%s<-%s\n' % (each, each_split_0_bin, each_split_1_bin))

    for each in overlap_list:
        output.write('%s\tN/A\n' % each)

    output.close()




add_direction(pwd_HGT_candidates_file, pwd_HGT_candidates_with_direction_file)





