import os
import shutil


usage = '''

module load python/3.5.2
python3 /srv/scratch/z5039045/Scripts/qsub_file_generator_COG_annotation.py

'''

##################################################### CONFIGURATION ####################################################

generate_qsub_file = 0
combine_results = 1
#sequence_file_list = '/Users/songweizhi/Desktop/faa_file_list.txt'
sequence_file_list = '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/faa_file_list.txt'


# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '11:59:00'
email = '244289990@qq.com'
modules_needed = ['python/3.5.2', 'perl/5.20.1', 'blast+/2.6.0']
qsub_file_folder = '/Users/songweizhi/Desktop/qsub_files_COG_annotation'
sequence_file_folder = '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/GoodBins_0.5_0.05_faa_files'


# parameters for combine annotation results
annotation_results_folder = '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/GoodBins_0.5_0.05_COG_annotation'
func_stats_file_folder =  '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/GoodBins_0.5_0.05_COG_annotation/0_func_stats_files'


################################################## generate_qsub_file ##################################################

if generate_qsub_file == 1:

    # Create qsub file folder
    if os.path.isdir(qsub_file_folder):
        shutil.rmtree(qsub_file_folder)
        if os.path.isdir(qsub_file_folder):
            shutil.rmtree(qsub_file_folder)
            if os.path.isdir(qsub_file_folder):
                shutil.rmtree(qsub_file_folder)
    os.system('mkdir %s' % qsub_file_folder)


    # Prepare header
    line_1 = '#!/bin/bash'
    line_2 = '#PBS -l nodes=%s:ppn=%s' % (str(nodes_number), str(ppn_number))
    line_3 = '#PBS -l vmem=%sgb' % str(memory)
    line_4 = '#PBS -l walltime=%s' % walltime_needed
    line_5 = '#PBS -j oe'
    line_6 = '#PBS -M %s' % email
    line_7 = '#PBS -m ae'
    header = '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % (line_1, line_2, line_3, line_4, line_5, line_6, line_7)


    # Prepare module lines
    module_lines = ''
    for module in modules_needed:
        module_lines += 'module load %s\n' % module


    # write to qsub files
    for each in open(sequence_file_list):
        current_qsub_file = '%s/qsub_COG_annotation_%s.sh' % (qsub_file_folder, each.strip())
        handle = open(current_qsub_file, 'w')
        handle.write('\n' + header + '\n\n')
        handle.write(module_lines)
        handle.write('\ncd %s\n' % sequence_file_folder)
        handle.write('python3 /srv/scratch/z5039045/Scripts/COG_wrapper.py -in %s -t P\n' % (each.strip()))
        handle.close()


############################################## combine annotation results ##############################################

if combine_results == 1:

    for each_file in open(sequence_file_list):
        each_file_no_ext = '.'.join(each_file.strip().split('.')[:-1])
        pwd_func_stats_file = '%s/%s_COG_results/func_stats.txt'   % (annotation_results_folder, each_file_no_ext)
        pwd_renamed_func_stats_file = '%s/%s_func_stats.txt'       % (func_stats_file_folder, each_file_no_ext)
        cmd_copy = 'cp %s %s' % (pwd_func_stats_file, pwd_renamed_func_stats_file)

        if os.path.isfile(pwd_func_stats_file) == 1:
            os.system(cmd_copy)
        else:
            print('Annotaion results not found, run again with: qsub qsub_COG_annotation_%s.sh' % each_file.strip())
