import os

out = open('/Users/songweizhi/Desktop/666/6_out.tab', 'w')
out_2 = open('/Users/songweizhi/Desktop/666/6_out_2.tab', 'w')

match_dict = {}
len_dict = {}

for each in open('/Users/songweizhi/Desktop/666/6.tab'):
    #print(each)
    each_split = each.strip().split('\t')
    iden = float(each_split[2])
    query = each_split[0].split('.')[0]
    subject = each_split[1]
    subject_new = subject.split('|')[3].split('.')[0]
    length = int(each_split[3])
    if (iden >= 100) and (length >= 30):
        #out.write(each)
        #print(each)
        out.write('%s\t%s\t%s\n' % (subject_new, query, length))
        match_dict[subject_new] = query
        len_dict[subject_new] = length

out.close()

os .system('cat %s | sort > %s' % ('/Users/songweizhi/Desktop/666/6_out.tab', '/Users/songweizhi/Desktop/666/6_out_sorted.tab'))

matched_num = 0
unmatched_num = 0
for each_g in open('/Users/songweizhi/Desktop/666/old_length.txt'):
    #print(each_g)
    id = each_g.strip().split('\t')[0]
    aa_len = each_g.strip().split('\t')[1]
    if id not in len_dict:
        print('%s\t%s' % ('-', each_g))
        out_2.write('%s\t%s' % ('-', each_g))
        unmatched_num += 1
    else:
        if len_dict[id] == int(aa_len):
            print('%s\t%s' % (match_dict[id], each_g))
            out_2.write('%s\t%s' % (match_dict[id], each_g))
            matched_num += 1
        else:
            print('%s\t%s' % ('-', each_g))
            out_2.write('%s\t%s' % ('-', each_g))
            unmatched_num += 1

out_2.close()
print(matched_num)
print(unmatched_num)






