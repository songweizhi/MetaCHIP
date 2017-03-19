
wd = '/Users/songweizhi/Desktop/new'
transfer_file = open('%s/donor2recip.txt' % wd)
predicted = open('%s/HGT_candidates.txt' % wd)


donor_list = []
for each in transfer_file:
    each_split = each.strip().split('\t')
    donor = each_split[0]
    if donor != 'donor_gene':
        donor_list.append(donor)

print('donor_list(%s):' % len(donor_list))
print('\t'.join(donor_list))

predicted_list = []
for each in predicted:
    each_split = each.strip().split('\t')
    each_split_sort = sorted(each_split)
    predicted_list.append(each_split_sort[0])

print('predicted_list(%s):' % len(predicted_list))
print('\t'.join(predicted_list))


validated = []
failed = []
for each in donor_list:
    if each in predicted_list:
        validated.append(each)
    if each not in predicted_list:
        failed.append(each)

print(validated)
print(failed)

print('validated(%s):' % len(validated))
print('\t'.join(validated))

print('failed(%s):' % len(failed))
print('\t'.join(failed))
