# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com
# t.thomas@unsw.edu.au

# Binning_refiner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Binning_refiner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import glob
import shutil
import argparse

"""

    It will:
    1. Create a new folder for each of your input bins and copy your bin into its corresponding folder
    2. Generate qsub file to run CheckM for each input bin
    3. Submit generated qsub files if specified

"""

#################################################### CONFIGURATION #####################################################

parser = argparse.ArgumentParser()

parser.add_argument('-1',
                    required=False,
                    help='first bin folder name')

parser.add_argument('-2',
                    required=False,
                    help='second bin folder name')

parser.add_argument('-3',
                    required=False,
                    help='third bin folder name')

parser.add_argument('-r',
                    required=False,
                    help='refined bin folder name')

parser.add_argument('-pbs',
                    action="store_true",
                    help='generate PBS job scripts')

parser.add_argument('-qsub',
                    action="store_true",
                    help='submit generated PBS job scripts')

parser.add_argument('-e',
                    required=False,
                    help='your email address',
                    metavar='(opt)')

parser.add_argument('-nodes',
                    required=False,
                    default=1,
                    type=int,
                    help='nodes number needed (default = 1)',
                    metavar='(opt)')

parser.add_argument('-ppn',
                    required=False,
                    default=1,
                    type=int,
                    help='ppn number needed (default = 1)',
                    metavar='(opt)')

parser.add_argument('-memory',
                    required=False,
                    default=120,
                    type=int,
                    help='memory needed (default = 120)',
                    metavar='(opt)')

parser.add_argument('-walltime',
                    required=False,
                    default='03:59:00',
                    help='walltime needed (default = 03:59:00)',
                    metavar='(opt)')

parser.add_argument('-python',
                    required=False,
                    default='python/2.7.8',
                    help='python version (default: python/2.7.8)',
                    metavar='(opt)')

parser.add_argument('-hmmer',
                    required=False,
                    default='hmmer/3.1b2',
                    help='hmmer version (default: hmmer/3.1b2)',
                    metavar='(opt)')

parser.add_argument('-pplacer',
                    required=False,
                    default='pplacer/1.1.alpha16',
                    help='pplacer version (default: pplacer/1.1.alpha16)',
                    metavar='(opt)')

parser.add_argument('-prodigal',
                    required=False,
                    default='prodigal/2.6.3',
                    help='prodigal version (default: prodigal/2.6.3)',
                    metavar='(opt)')


args = vars(parser.parse_args())

input_bin_folder_list = []

if args['1'] != None:
    input_bin_folder_1 = args['1']
    if input_bin_folder_1[-1] == '/':
        input_bin_folder_1 = input_bin_folder_1[:-1]
    input_bin_folder_list.append(input_bin_folder_1)

if args['2'] != None:
    input_bin_folder_2 = args['2']
    if input_bin_folder_2[-1] == '/':
        input_bin_folder_2 = input_bin_folder_2[:-1]
    input_bin_folder_list.append(input_bin_folder_2)

if args['3'] != None:
    input_bin_folder_3 = args['3']
    if input_bin_folder_3[-1] == '/':
        input_bin_folder_3 = input_bin_folder_3[:-1]
    input_bin_folder_list.append(input_bin_folder_3)

if args['r'] != None:
    input_refined_bin_folder = args['r']
    if input_refined_bin_folder[-1] == '/':
        input_refined_bin_folder = input_refined_bin_folder[:-1]
    input_bin_folder_list.append(input_refined_bin_folder)

print('Specified %s input bin sets: %s' % (len(input_bin_folder_list), ' '.join(input_bin_folder_list)))

nodes_number = args['nodes']
ppn_number = args['ppn']
memory = args['memory']
walltime_needed = args['walltime']
email = args['e']
modules_needed = [args['python'], args['hmmer'], args['pplacer'], args['prodigal']]

########################################################################################################################

# define folder/file name
wd = os.getcwd()
checkm_wd = 'checkm_wd_universal'

# prepare qsub file header and module lines
header = ''
module_lines = ''
if args['pbs'] == True:
    # prepare qsub file header
    line_1 = '#!/bin/bash\n'
    line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
    line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
    line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
    line_5 = '#PBS -j oe\n'
    line_6 = ''
    if args['e'] != None:
        line_6 = '#PBS -M ' + email + '\n'
    line_7 = '#PBS -m ae\n'
    line_8 = 'cd $PBS_O_WORDIR\n'
    #line_8 = 'cd %s\n' % pwd_qsub_files_folder
    header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8

    # Prepare module lines
    module_lines = ''
    for module in modules_needed:
        module_lines += 'module load ' + module + '\n'

for each_bin_folder in input_bin_folder_list:
    pwd_checkm_wd = '%s/%s/%s' % (wd, each_bin_folder, checkm_wd)

    choice = ''
    if os.path.isdir(pwd_checkm_wd):
        choice = str(input('CheckM working directory detected, Press "S/s" to skip this step, Press any other key to overwrite it.\nYour choice: '))

    if choice in ['S', 's']:
        pass
    else:
        # get CheckM working directory folder
        if os.path.isdir(pwd_checkm_wd):
            shutil.rmtree(pwd_checkm_wd)
            os.mkdir(pwd_checkm_wd)
        else:
            os.mkdir(pwd_checkm_wd)

        # get bin name list
        bin_files = '%s/%s/*.fa*' % (wd, each_bin_folder)
        bin_list = [os.path.basename(file_name) for file_name in glob.glob(bin_files)]

        if len(bin_list) == 0:
            print('No input bin detected from %s, please double-check.' % pwd_checkm_wd)
            exit()

        bin_file_ext_list = []
        for bin in bin_list:
            name, ext = os.path.splitext(bin)
            bin_file_ext_list.append(ext[1:])

        # uniq bin_file_ext_list
        bin_file_ext_list_uniq = []
        for each in bin_file_ext_list:
            if each not in bin_file_ext_list_uniq:
                bin_file_ext_list_uniq.append(each)
            else:
                pass

        # check whether bins in the same folder have same extension, exit if not
        if len(bin_file_ext_list_uniq) > 1:
            print('Different file extensions were foud from %s bins, please use same extension (fa, fas or fasta) '
                  'for all bins in the same folder.' % each_bin_folder)
            exit()
        else:
            pass

        # get bin file extension
        bin_file_extension = bin_file_ext_list_uniq[0]

        # get qsub file and submit it
        qsub_files_folder = 'qsub_files_for_checkm'
        pwd_qsub_files_folder = '%s/%s/%s/%s' % (wd, each_bin_folder, checkm_wd, qsub_files_folder)
        if args['pbs'] == True:
            os.mkdir(pwd_qsub_files_folder)

        if '/' in each_bin_folder:
            each_bin_folder_split = each_bin_folder.split('/')
            checkm_cmds_file_name = 'Commands_for_CheckM_%s.txt' % each_bin_folder[-1]
        else:
            checkm_cmds_file_name = 'Commands_for_CheckM_%s.txt' % each_bin_folder

        pwd_checkm_cmds_file_name = '%s/%s' % (pwd_checkm_wd, checkm_cmds_file_name)
        #checkm_cmds_file = open(pwd_checkm_cmds_file_name, 'w')

        for each_bin in bin_list:
            print(each_bin)
            bin_name = each_bin[:-(len(bin_file_extension) + 1)]

            # create a folder for current bin
            checkm_wd_for_current_bin= '%s/%s/%s/%s' % (wd, each_bin_folder, checkm_wd, bin_name)
            os.mkdir(checkm_wd_for_current_bin)
            pwd_bin = '%s/%s/%s' % (wd, each_bin_folder, each_bin)
            os.system('cp %s %s' % (pwd_bin, checkm_wd_for_current_bin))

            if args['pbs'] == True:
                qsub_file = '%s.sh' % bin_name
                pwd_qsub_file = '%s/%s' % (pwd_qsub_files_folder, qsub_file)
                qsub_file_handle = open(pwd_qsub_file, 'w')
                qsub_file_handle.write('%s\n%s' % (header, module_lines))

                cmd_analyze = 'checkm analyze /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm %s %s -x %s\n' % (checkm_wd_for_current_bin,checkm_wd_for_current_bin, bin_file_extension)
                cmd_qa = 'checkm qa /srv/scratch/z5039045/Softwares/MetaCHIP/phylo.hmm %s -f %s/%s.txt\n' % (checkm_wd_for_current_bin, pwd_checkm_wd, bin_name)

                cmds = cmd_analyze + cmd_qa

                #cmds =          'checkm lineage_wf -x %s -t %s %s %s -f %s/%s.txt\n'   % (bin_file_extension, ppn_number, checkm_wd_for_current_bin, checkm_wd_for_current_bin, pwd_checkm_wd, bin_name)
                #cmds_non_qsub = 'checkm lineage_wf -x %s -t %s %s %s -f %s/%s.txt &\n' % (bin_file_extension, ppn_number, checkm_wd_for_current_bin, checkm_wd_for_current_bin, pwd_checkm_wd, bin_name)

                qsub_file_handle.write(cmds)
                #checkm_cmds_file.write(cmds_non_qsub)
                qsub_file_handle.close()

                if args['qsub'] == True:
                    os.chdir(pwd_qsub_files_folder)
                    os.system('qsub %s' % pwd_qsub_file)
                    os.chdir(wd)

        #checkm_cmds_file.close()
