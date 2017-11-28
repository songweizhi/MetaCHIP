import os

mutation_group = 'm30'

output_handle = open('/Users/songweizhi/Desktop/ctg_assign/%s_ctg_assignment_id.tab' % mutation_group, 'w')
for each in open('/Users/songweizhi/Desktop/ctg_assign/%s_ctg_assignment.tab' % mutation_group):
    #print(each)
    each_split = each.strip().split('\t')
    iden = float(each_split[2])
    align_len = int(each_split[3])
    query_len = int(each_split[12])
    subject_len = int(each_split[13])
    query_cov = align_len/query_len*100
    if (iden >= 99) and query_cov >= 99:
        print('%s\t%s' % (each_split[0], each_split[1]))
        output_handle.write('%s\t%s\n' % (each_split[0], each_split[1]))
output_handle.close()

os.system('cat /Users/songweizhi/Desktop/ctg_assign/%s_ctg_assignment_id.tab | sort | uniq > /Users/songweizhi/Desktop/ctg_assign/%s_ctg_assignment_id_uniq.tab' % (mutation_group, mutation_group))
os.remove('/Users/songweizhi/Desktop/ctg_assign/%s_ctg_assignment_id.tab' % mutation_group)
