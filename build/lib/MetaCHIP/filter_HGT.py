import os
import shutil


filter_HGT_usage = '''
====================================== filter_HGT example commands ======================================

# get HGTs detected at at least TWO levels
MetaCHIP filter_HGT -i NorthSea_pcofg_detected_HGTs.txt -n 2

# get HGTs detected at at least THREE levels and copy their flanking region plots into a new folder
MetaCHIP filter_HGT -i NorthSea_pcofg_detected_HGTs.txt -n 3 -plot NorthSea_pcofg_Flanking_region_plots

=========================================================================================================
'''


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def filter_HGT(args):

    file_in =           args['i']
    n =                 args['n']
    flk_plot_folder =   args['plot']

    file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
    file_out = '%s/%s_min_level_num_%s%s' % (file_in_path, file_in_basename, n, file_in_extension)

    flk_plot_folder_qualified = None
    if flk_plot_folder is not None:

        if flk_plot_folder[-1] == '/':
            flk_plot_folder = flk_plot_folder[:-1]

        flk_plot_folder_path = '.'
        flk_plot_folder_name = flk_plot_folder
        if '/' in flk_plot_folder:
            flk_plot_folder_path = '/'.join(flk_plot_folder.split('/')[:-1])
            flk_plot_folder_name = flk_plot_folder.split('/')[-1]

        flk_plot_folder_qualified = '%s/%s_min_level_num_%s' % (flk_plot_folder_path, flk_plot_folder_name, n)
        force_create_folder(flk_plot_folder_qualified)

    file_out_handle = open(file_out, 'w')
    file_out_handle.write(open(file_in).readline())
    if 'occurence' not in open(file_in).readline().strip():
        print('Not multiple level predictions, filter_HGT exited')
        exit()
    else:
        for predicted_hgt in open(file_in):

            if not predicted_hgt.startswith('Gene_1	Gene_2	Identity'):
                predicted_hgt_split = predicted_hgt.strip().split('\t')
                hgt_occurence = predicted_hgt_split[3]
                hgt_occurence_1_num = hgt_occurence.count('1')

                if hgt_occurence_1_num >= n:
                    file_out_handle.write(predicted_hgt)

                    if flk_plot_folder is not None:
                        pwd_hgt_plot = '%s/%s___%s.SVG' % (flk_plot_folder, predicted_hgt_split[0], predicted_hgt_split[1])
                        cp_cmd = 'cp %s %s/' % (pwd_hgt_plot, flk_plot_folder_qualified)
                        os.system(cp_cmd)

    file_out_handle.close()
