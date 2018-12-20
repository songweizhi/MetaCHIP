
pwd_phylo_hmm_stat_txt = '/Users/songweizhi/Desktop/phylo.hmm.stat.txt'

pwd_hmm_id_file = ''
pwd_hmm_id_file_handle = open(pwd_hmm_id_file, 'w')
for each_profile in open(pwd_phylo_hmm_stat_txt):
    if not each_profile.startswith('#'):
        each_profile_split = each_profile.strip().split(' ')
        if each_profile_split != ['']:
            each_profile_split_no_space = []
            for each_element in each_profile_split:
                if each_element != '':
                    each_profile_split_no_space.append(each_element)
            hmm_profile_id = each_profile_split_no_space[2]
            pwd_hmm_id_file_handle.write('%s\n' % hmm_profile_id)
pwd_hmm_id_file_handle.close()

