

pwd_BM_HGTs = '/Users/songweizhi/Desktop/MetaCHIP_test/NorthSea_MetaCHIP_wd/NorthSea_c5_ip90_al200bp_c70_ei95bp_f10kbp/NorthSea_c5_HGTs_BM.txt'
pwd_blast_subjects_in_one_line = '/Users/songweizhi/Desktop/MetaCHIP_test/NorthSea_MetaCHIP_wd/NorthSea_c5_ip90_al200bp_c70_ei95bp_f10kbp/NorthSea_c5_subjects_in_one_line.txt'
pwd_query_to_subjects_file = '/Users/songweizhi/Desktop/MetaCHIP_test/NorthSea_MetaCHIP_wd/NorthSea_c5_ip90_al200bp_c70_ei95bp_f10kbp/NorthSea_c5_HGT_subjects.txt'

def export_HGT_query_to_subjects(pwd_BM_HGTs, pwd_blast_subjects_in_one_line, pwd_query_to_subjects_file):

    HGT_candidates = set()
    for HGT_pair in open(pwd_BM_HGTs):
        HGT_pair_split = HGT_pair.strip().split('\t')
        gene_1 = HGT_pair_split[0]
        gene_2 = HGT_pair_split[1]
        HGT_candidates.add(gene_1)
        HGT_candidates.add(gene_2)

    query_subjects_dict = {}
    for each_gene in open(pwd_blast_subjects_in_one_line):
        each_gene_split = each_gene.strip().split('\t')
        query = each_gene_split[0].split('|')[1]
        subjects = [i.split('|')[1] for i in each_gene_split[1:]]
        if query in HGT_candidates:
            query_subjects_dict[query] = subjects

    pwd_query_to_subjects_file_handle = open(pwd_query_to_subjects_file, 'w')
    for each in query_subjects_dict:
        for_out = '%s\t%s\n' % (each, ','.join(query_subjects_dict[each]))
        pwd_query_to_subjects_file_handle.write(for_out)
    pwd_query_to_subjects_file_handle.close()


export_HGT_query_to_subjects(pwd_BM_HGTs, pwd_blast_subjects_in_one_line, pwd_query_to_subjects_file)
