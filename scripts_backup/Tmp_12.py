import os
import shutil




##################################################### CONFIGURATION ####################################################

# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 3
memory = 20
walltime_needed = '11:59:00'
email = 'weizhi.song@student.unsw.edu.au'
modules_needed = ['blast+/2.6.0']
job_script_folder = '/Users/songweizhi/Desktop/qsub_files_blastn'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_Hung'


################################################## generate_qsub_file ##################################################


def create_job_script(job_script_folder, job_script_file_name, node_num, ppn_num, memory, walltime, email, modules_list, cmd, qsub):

    # create job script folder
    if os.path.isdir(job_script_folder):
        shutil.rmtree(job_script_folder)
        if os.path.isdir(job_script_folder):
            shutil.rmtree(job_script_folder)
            if os.path.isdir(job_script_folder):
                shutil.rmtree(job_script_folder)
    os.system('mkdir %s' % job_script_folder)


    # Prepare header
    line_1 = '#!/bin/bash'
    line_2 = '#PBS -l nodes=%s:ppn=%s' % (str(node_num), str(ppn_num))
    line_3 = '#PBS -l vmem=%sgb' % str(memory)
    line_4 = '#PBS -l walltime=%s' % walltime
    line_5 = '#PBS -j oe'
    line_6 = '#PBS -M %s' % email
    line_7 = '#PBS -m ae'
    header = '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % (line_1, line_2, line_3, line_4, line_5, line_6, line_7)


    # Prepare module lines
    module_lines = ''
    for module in modules_list:
        module_lines += 'module load %s\n' % module


    # write to qsub files
    output_file_handle = open('%s/%s' % (job_script_folder, job_script_file_name), 'w')
    output_file_handle.write(header)
    output_file_handle.write(module_lines)
    output_file_handle.write('cd %s\n' % wd_on_katana)
    output_file_handle.write(cmd)
    output_file_handle.close()
