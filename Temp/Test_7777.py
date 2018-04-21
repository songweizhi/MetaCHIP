
total = 0
num_B = []
num_D = []
num_H = []
num_L = []
num_M = []
for each in open('/Users/songweizhi/Desktop/HGT_candidates.txt'):
    each_split = each.strip().split('\t')
    gene_1_group = each_split[2]
    gene_2_group = each_split[3]
    if (gene_1_group[0] == 'B') and (gene_2_group[0] == 'B'):
        num_B.append('%s_%s' % (gene_1_group, gene_2_group))
    if (gene_1_group[0] == 'D') and (gene_2_group[0] == 'D'):
        num_D.append('%s_%s' % (gene_1_group, gene_2_group))
    if (gene_1_group[0] == 'H') and (gene_2_group[0] == 'H'):
        num_H.append('%s_%s' % (gene_1_group, gene_2_group))
    if (gene_1_group[0] == 'L') and (gene_2_group[0] == 'L'):
        num_L.append('%s_%s' % (gene_1_group, gene_2_group))
    if (gene_1_group[0] == 'M') and (gene_2_group[0] == 'M'):
        num_M.append('%s_%s' % (gene_1_group, gene_2_group))

    total += 1

print(total)
print(num_B)
print(num_D)
print(num_H)
print(num_L)
print(num_M)



print('B:\t%s\t%s' % (len(num_B), len(num_B)/total*100))
print('D:\t%s\t%s' % (len(num_D), len(num_D)/total*100))
print('H:\t%s\t%s' % (len(num_H), len(num_H)/total*100))
print('L:\t%s\t%s' % (len(num_L), len(num_L)/total*100))
print('M:\t%s\t%s' % (len(num_M), len(num_M)/total*100))


print('%s(%s)' % (len(num_B), float("{0:.2f}".format(len(num_B)/total*100))))
print('%s(%s)' % (len(num_D), float("{0:.2f}".format(len(num_D)/total*100))))
print('%s(%s)' % (len(num_H), float("{0:.2f}".format(len(num_H)/total*100))))
print('%s(%s)' % (len(num_L), float("{0:.2f}".format(len(num_L)/total*100))))
print('%s(%s)' % (len(num_M), float("{0:.2f}".format(len(num_M)/total*100))))

total_within_group_HGTs = len(num_B) + len(num_D) + len(num_H) + len(num_L) + len(num_M)
total_within_group_HGTs_percent = float("{0:.2f}".format(total_within_group_HGTs/total*100))


print(total_within_group_HGTs)
print(total_within_group_HGTs_percent)