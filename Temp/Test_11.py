# import os
# import matplotlib.pyplot as plt
# import pylab
#
# blast_results = '/Users/songweizhi/Desktop/all_vs_all_ref.tab'
# new_1_file = '/Users/songweizhi/Desktop/identities_1.txt'
# #new_2_file = '/Users/songweizhi/Desktop/identities_2.txt'
#
# new_1_handle = open(new_1_file, 'w')
# identities = []
# for each in open(blast_results):
#     #print(each)
#     each_split = each.strip().split('\t')
#     identity = float(each_split[2])
#     length = int(each_split[3])
#     if identity >= 95 and 100000 <= length < 3758219:
# #    if identity >= 95 and 1000 <= length < 517008:
#         new_1_handle.write('%s\n' % length)
#         identities.append(length)
#
# identities_sorted = sorted(identities)
# #os.system('cat %s | sort > %s' % (new_1_file, new_2_file))
# print(identities)
# print(identities_sorted)
# print(len(identities_sorted))
#
#
# plt.hist(identities_sorted)
# plt.show()
n = 0
for each in open('/Users/songweizhi/Desktop/double_check/0_vs_5_ffn.tab'):
    #print(each)
    each_split = each.strip().split('\t')
    identity = float(each_split[2])
    query_len = each_split[12]
    subject_len = each_split[13]
    align_len = each_split[3]
    if each_split[0] == each_split[1]:
        print(each)
    # if identity == 100 and query_len == align_len:
    #     print(each)
        n += 1

print(n)
