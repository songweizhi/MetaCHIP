import os
import glob


def rm_folder_file(target_re):
    target_list = glob.glob(target_re)

    for target in target_list:

        if os.path.isdir(target) is True:
            os.system('rm -r %s' % target)

        elif os.path.isfile(target) is True:
            os.system('rm %s' % target)


act_file_re =   'Subsample_2000_f298_Flanking_region_plots_tmp/*___*/*'
act_folder_re = 'Subsample_2000_f298_Flanking_region_plots_tmp/*___*'

rm_folder_file(act_file_re)
rm_folder_file(act_folder_re)
