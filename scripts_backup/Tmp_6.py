


selected_gene_ortholog_dict = {}
filtered = '/Users/songweizhi/Desktop/kerrin/kerrin_selection_qualified.tab'
filtered_handle = open(filtered, 'w')

for match in open('/Users/songweizhi/Desktop/kerrin/kerrin_selection.tab'):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    query_genome = '_'.join(query.split('_')[:-1])
    subject_genome = '_'.join(subject.split('_')[:-1])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))

    if (query != subject) and (coverage_q >= 80) and (coverage_s >= 80):

        #print(match.strip())
        filtered_handle.write(match)

        if query not in selected_gene_ortholog_dict:
            selected_gene_ortholog_dict[query] = [set(), set()]

        if subject_genome.startswith('A'):
            selected_gene_ortholog_dict[query][0].add(subject_genome)

        if subject_genome.startswith('B'):
            selected_gene_ortholog_dict[query][1].add(subject_genome)

filtered_handle.close()

print(len(selected_gene_ortholog_dict))


for each in selected_gene_ortholog_dict:
    print('%s\t%s\t%s' % (each, len(selected_gene_ortholog_dict[each][0]), len(selected_gene_ortholog_dict[each][1])))



