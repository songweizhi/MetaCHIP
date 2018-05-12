import os


def get_hits_group(input_file_name, output_file_name):
    matches_2 = open(input_file_name)
    output_2_file = open(output_file_name, 'w')
    current_gene = ''
    group_member = []
    for match in matches_2:
        match_split = match.strip().split('\t')
        query_2 = match_split[0]
        target_2 = match_split[1]
        if current_gene == '':
            current_gene = query_2
            group_member.append(target_2)
        else:
            if query_2 == current_gene:
                if target_2 not in group_member:
                    group_member.append(target_2)
                else:
                    pass
            else:
                output_2_file.write('%s\t%s' % (current_gene, '\t'.join(group_member)) + '\n')
                current_gene = query_2
                group_member = []
                group_member.append(target_2)
    output_2_file.write('%s\t%s' % (current_gene, '\t'.join(group_member)) + '\n')
    output_2_file.close()


wd = '/Users/songweizhi/Desktop/get_core'
pwd_matches = '%s/%s' % (wd, 'all_vs_all_ffn.tab')
pwd_qualified_matches = '%s/%s' % (wd, 'qualified_matches.txt')
pwd_qualified_matches_sorted = '%s/%s' % (wd, 'qualified_matches_sorted.txt')
pwd_qualified_matches_sorted_one_line = '%s/%s' % (wd, 'qualified_matches_sorted_one_line.txt')
pwd_qualified_matches_sorted_one_line_reorded = '%s/%s' % (wd, 'qualified_matches_sorted_one_line_reorded.txt')
pwd_qualified_matches_sorted_one_line_reorded_sorted_uniq = '%s/%s' % (wd, 'qualified_matches_sorted_one_line_reorded_sorted.txt')
pwd_qualified_core_genes = '%s/%s' % (wd, 'qualified_core_genes.txt')

matches = open(pwd_matches)
qualified_matches = open(pwd_qualified_matches, 'w')

for match in matches:
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    query_split = query.split('_')
    subject_split = subject.split('_')
    query_bin_name = '_'.join(query_split[:-1])
    subject_bin_name = '_'.join(subject_split[:-1])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))
    if (coverage_q >= 80) and (coverage_s >= 80) and (identity >= 70) and (align_len >= 500):
        qualified_matches.write('%s\t%s\n' % (query, subject))
qualified_matches.close()

os.system('cat %s | sort > %s' % (pwd_qualified_matches, pwd_qualified_matches_sorted))

get_hits_group(pwd_qualified_matches_sorted, pwd_qualified_matches_sorted_one_line)


pwd_qualified_matches_sorted_one_line_handle = open(pwd_qualified_matches_sorted_one_line)
pwd_qualified_matches_sorted_one_line_reorded_handle = open(pwd_qualified_matches_sorted_one_line_reorded, 'w')
for each in pwd_qualified_matches_sorted_one_line_handle:
    each_split = each.strip().split('\t')
    each_split_uniq = []
    for each in each_split:
        if each not in each_split_uniq:
            each_split_uniq.append(each)

    each_split_uniq_sorted = sorted(each_split_uniq)
    if len(each_split_uniq_sorted) > 1:
        pwd_qualified_matches_sorted_one_line_reorded_handle.write('\t'.join(each_split_uniq_sorted) + '\n')
pwd_qualified_matches_sorted_one_line_reorded_handle.close()

os.system('cat %s | sort | uniq > %s' % (pwd_qualified_matches_sorted_one_line_reorded, pwd_qualified_matches_sorted_one_line_reorded_sorted_uniq))

pwd_qualified_matches_sorted_one_line_reorded_sorted_uniq_handle = open(pwd_qualified_matches_sorted_one_line_reorded_sorted_uniq)
pwd_qualified_core_genes_handle = open(pwd_qualified_core_genes, 'w')

exported_group = []
for each in pwd_qualified_matches_sorted_one_line_reorded_sorted_uniq_handle:
    each_split = each.strip().split('\t')
    each_split_A = []
    each_split_B = []
    for each_2 in each_split:
        if each_2[0] == 'A':
            each_split_A.append(each_2)
        if each_2[0] == 'B':
            each_split_B.append(each_2)
    if (len(each_split_A) >= 3) and (len(each_split_B) >= 3) and (each_split[0] not in exported_group):
        pwd_qualified_core_genes_handle.write(each)
        exported_group.append(each_split[0])
pwd_qualified_core_genes_handle.close()


########################################################################################################################

transfers = open('%s/donor2recip.txt' % wd)

selected = []
chongfude = []
for each in transfers:
    each_split_n = each.strip().split('\t')
    print(each_split_n[0])
    if each_split_n[0] not in selected:
        selected.append(each_split_n[0])
    else:
        chongfude.append(each_split_n[0])

print(len(selected))
print(chongfude)

cores = open(pwd_qualified_core_genes)
for each_core in cores:
    each_core_split = each_core.strip().split('\t')
    existed = 0
    for each_gene in each_core_split:
        if each_gene in selected:
            existed = 1
    if existed == 0:
        print(each_core.strip())
