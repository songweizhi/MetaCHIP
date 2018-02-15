import os

# script to get the average coverage for each bin from the checkm coverage output file
# the first line need to be removed from the output file
# currently only support for three replicates, that is three bam files

os.chdir('/Users/songweizhi/Desktop/get_average_coverage')


out = open('coverage_each_bin.tsv', 'w')
for each in open('coverage.tsv'):
    each_split = each.strip().split('\t')
    bin_id = each_split[1]
    sample_1 = each_split[4]
    sample_2 = each_split[7]
    sample_3 = each_split[10]
    for_write = '%s\t%s\t%s\t%s\n' % (bin_id, sample_1, sample_2, sample_3)
    if bin_id != 'unbinned':
        out.write(for_write)
out.close()


os.system('cat coverage_each_bin.tsv | sort > coverage_each_bin_sorted.tsv')


out2 = open('coverage_each_bin_average', 'w')
current_bin = ''
ctg_number = 0
sample1_sum = 0
sample2_sum = 0
sample3_sum = 0
for each_ctg in open('coverage_each_bin_sorted.tsv'):
    each_ctg_split = each_ctg.strip().split('\t')
    bin = each_ctg_split[0]
    sample1 = float(each_ctg_split[1])
    sample2 = float(each_ctg_split[2])
    sample3 = float(each_ctg_split[3])
    if current_bin == '':
        current_bin = bin
        ctg_number += 1
        sample1_sum += sample1
        sample2_sum += sample2
        sample3_sum += sample3
    elif current_bin == bin:
        ctg_number += 1
        sample1_sum += sample1
        sample2_sum += sample2
        sample3_sum += sample3
    elif current_bin != bin:
        out2.write('%s\t%s\t%s\t%s\n' % (current_bin,
                                         float("{0:.3f}".format(sample1_sum/ctg_number)),
                                         float("{0:.3f}".format(sample2_sum/ctg_number)),
                                         float("{0:.3f}".format(sample3_sum/ctg_number))))
        current_bin = bin
        ctg_number = 1
        sample1_sum = sample1
        sample2_sum = sample2
        sample3_sum = sample3
out2.write('%s\t%s\t%s\t%s\n' % (current_bin,
                                 float("{0:.3f}".format(sample1_sum/ctg_number)),
                                 float("{0:.3f}".format(sample2_sum/ctg_number)),
                                 float("{0:.3f}".format(sample3_sum/ctg_number))))
out2.close()

