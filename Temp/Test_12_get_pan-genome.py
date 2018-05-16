import os

input = '/Users/songweizhi/Desktop/D6Cluster_f70_1taxa_algCOG_e0_E1e-05_C70_S70_.cluster_list'
output = '/Users/songweizhi/Desktop/sumary.txt'

def get_element_number(list, wanted):
    num = 0
    for each in list:
        if each == wanted:
            num += 1
    return num

current_cluster = ''
members = []
output_handle = open(output, 'w')
output_handle.write('Cluster\tB06_Cluster.2\tD6_Cluster.6\tH8_Cluster.1\tH13_metabat_bin.6\n')
for each in open(input):
    each = each.strip()
    if each.startswith('cluster'):
        cluster = each.split(' ')[1]
        if current_cluster == '':
            current_cluster = cluster
        if current_cluster != '' and cluster != current_cluster:
            output_handle.write('%s\t%s\t%s\t%s\t%s\n' % (current_cluster, get_element_number(members, 'B06_Cluster.2'), get_element_number(members, 'D6_Cluster.6'), get_element_number(members, 'H8_Cluster.1'), get_element_number(members, 'H13_metabat_bin.6')))
            current_cluster = cluster
            members = []
    elif each.startswith(':'):
        current_member = each.split(' ')[1]
        genome_name, file_ext = os.path.splitext(current_member)
        members.append(genome_name)

output_handle.write('%s\t%s\t%s\t%s\t%s\n' % (
current_cluster, get_element_number(members, 'B06_Cluster.2'), get_element_number(members, 'D6_Cluster.6'),
get_element_number(members, 'H8_Cluster.1'), get_element_number(members, 'H13_metabat_bin.6')))
output_handle.close()

sum_1 = 0
sum_2 = 0
sum_3 = 0
sum_4 = 0
sum_12 = 0
sum_13 = 0
sum_14 = 0
sum_23 = 0
sum_24 = 0
sum_34 = 0
sum_123 = 0
sum_124 = 0
sum_134 = 0
sum_234 = 0
sum_1234 = 0
for each in open(output):
    if not each.startswith('Cluster'):
        each_split = each.strip().split('\t')
        num_1 = int(each_split[1])
        num_2 = int(each_split[2])
        num_3 = int(each_split[3])
        num_4 = int(each_split[4])
        if num_1 > 0:
            sum_1 += 1
        if num_2 > 0:
            sum_2 += 1
        if num_3 > 0:
            sum_3 += 1
        if num_4 > 0:
            sum_4 += 1
        if num_1 > 0 and num_2 > 0:
            sum_12 += 1
        if num_1 > 0 and num_3 > 0:
            sum_13 += 1
        if num_1 > 0 and num_4 > 0:
            sum_14 += 1
        if num_2 > 0 and num_3 > 0:
            sum_23 += 1
        if num_2 > 0 and num_4 > 0:
            sum_24 += 1
        if num_3 > 0 and num_4 > 0:
            sum_34 += 1
        if num_1 > 0 and num_2 > 0 and num_3 > 0:
            sum_123 += 1
        if num_1 > 0 and num_2 > 0 and num_4 > 0:
            sum_124 += 1
        if num_1 > 0 and num_3 > 0 and num_4 > 0:
            sum_134 += 1
        if num_2 > 0 and num_3 > 0 and num_4 > 0:
            sum_234 += 1
        if num_1 > 0 and num_2 > 0 and num_3 > 0 and num_4 > 0:
            sum_1234 += 1

print('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' % (sum_1,sum_2,sum_3,sum_4,sum_12,sum_13,sum_14,sum_23,sum_24,sum_34,sum_123,sum_124,sum_134,sum_234,sum_1234))
