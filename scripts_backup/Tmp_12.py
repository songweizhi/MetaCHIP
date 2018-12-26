out_handle = open('/Users/songweizhi/Desktop/combined_quality_filtered.txt', 'w')

n = 0
for each in open('/Users/songweizhi/Desktop/combined_quality.txt'):
    if not each.startswith('-'):
        if not each.startswith('  Bin Id'):
            #print(each)
            each_split = each.strip().split(' ')
            #print(each_split)

            each_split_no_space = []
            for each_element in each_split:

                if each_element != '':
                    each_split_no_space.append(each_element)
            #print(each_split_no_space)

            bin_id = each_split_no_space[0]
            completeness = float(each_split_no_space[12])
            contamination = float(each_split_no_space[13])
            heterogeneity = float(each_split_no_space[14])

            #print('%s\t%s\t%s\t%s' % (bin_id, completeness, contamination, heterogeneity))
            if (completeness >= 40) and (contamination <= 5):
                print('cp BetterBins_noNs_ctg2kbp_bin200kbp/%s.fa MetaBAT_bins_good/ &' % bin_id)
                n += 1

            out_handle.write(each)
out_handle.close()

print(n)









