
metachip_output = '/Users/songweizhi/Desktop/HGT_candidates_ET_validated.txt'

list_0 = 0
list_5 = 0
list_10 = 0
list_15 = 0
list_20 = 0
list_25 = 0
list_30 = 0
for each in open(metachip_output):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        identity = float(each_split[4])
        if 100 >= identity > 97.5:
            list_0 += 1
        elif 97.5 >= identity > 92.5:
            list_5 += 1
        elif 92.5 >= identity > 87.5:
            list_10 += 1
        elif 87.5 >= identity > 82.5:
            list_15 += 1
        elif 82.5 >= identity > 77.5:
            list_20 += 1
        elif 77.5 >= identity > 72.5:
            list_25 += 1
        else:
            list_30 += 1

print('Identity distribution:')
print('Mutation\t0\t5\t10\t15\t20\t25\t30\nHGT_number\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (list_0, list_5, list_10, list_15, list_20, list_25, list_30))








