import os
import glob


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


##################################################### CONFIGURATION ####################################################

# cd to wd_on_katana
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/between_genus/MetaCHIP_wd'
#wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_rebuttal/between_class/MetaCHIP_wd'

mutation_levels = [0, 5, 10, 15, 20, 25, 30]

################################################## generate_qsub_file ##################################################

# write to qsub files
bootstrap_num = 1
while bootstrap_num <= 10:
    for mutation_level in mutation_levels:
        get_homo_wd = 'bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_gbk_files_homologues' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        homo_folder_re = '%s/*taxa_*_S*_' % get_homo_wd
        homologues_folder = '%s/%s' % (get_homo_wd, [os.path.basename(file_name) for file_name in glob.glob(homo_folder_re)][0])
        homologues_folder_new = 'bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_homologues' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        rename_cmd = 'mv %s %s' % (homologues_folder, homologues_folder_new)

        print(rename_cmd)
        print('rm -r %s' % get_homo_wd)

        # os.system(rename_cmd)
        # os.system('rm -r %s' % get_homo_wd)

    bootstrap_num += 1
