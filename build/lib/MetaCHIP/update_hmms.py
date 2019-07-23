import os
import argparse
from time import sleep
from datetime import datetime


update_hmms_usage = '''
====================================== update_hmms example commands ====================================

# Download Pfam DB file (e.g. v32.0)
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# Download TIGRFAM DB folder (e.g. v14.0)
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz
tar -xzvf TIGRFAMs_14.0_HMM.tar.gz

MetaCHIP update_hmms -hmm MetaCHIP_phylo.hmm -p_db Pfam-A.hmm -t_db TIGRFAMs_14.0_HMM

========================================================================================================
'''

def update_hmms(args):

    current_43_hmm_file =        args['hmm']
    downloaded_pfam_db =         args['p_db']
    TIGRFAMs_profiles_folder =   args['t_db']

    #################################### define file name ####################################

    updated_43_profiles =       'MetaCHIP_phylo_updated.hmm'
    update_log_file =           'update_log.txt'
    current_43_hmm_file_stat =  'MetaCHIP_phylo.hmm.stat'
    downloaded_pfam_db_stats =  'Pfam-A.hmm.stat'
    Pfam_41_ids_updated =       'Pfam_41_ids_with_updated_version.txt'
    Pfam_profiles_updated =     'Pfam_profiles_updated.hmm'

    time_format = '[%Y-%m-%d %H:%M:%S]'

    ########################################## main ##########################################

    # get MetaCHIP_phylo.hmm stat
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Get summary statistics for %s' % current_43_hmm_file))
    hmmstat_cmd_MetaCHIP = 'hmmstat %s > %s' % (current_43_hmm_file, current_43_hmm_file_stat)
    os.system(hmmstat_cmd_MetaCHIP)

    # read in 41 Pfam ids
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Parse summary statistics for %s' % current_43_hmm_file))
    Pfam_41_id_list_with_version = []
    Pfam_41_id_list_no_version = []
    TIGRFAMs_2_id_list = []
    id_to_version_dict_old = {}
    for each_hmm in open(current_43_hmm_file_stat):

        if not each_hmm.startswith('#') and (each_hmm.strip() != ''):
            each_hmm_split = each_hmm.strip().split(' ')

            # remove '' from each_hmm_split
            each_hmm_split_no_space = []
            for i in each_hmm_split:
                if i != '':
                    each_hmm_split_no_space.append(i)

            hmm_id = each_hmm_split_no_space[2]
            if hmm_id.startswith('PF'):
                Pfam_41_id_list_with_version.append(hmm_id)
                Pfam_41_id_list_no_version.append(hmm_id.split('.')[0])
                id_to_version_dict_old[hmm_id.split('.')[0]] = hmm_id
            else:
                TIGRFAMs_2_id_list.append(hmm_id)


    # get all hmm profile id from downloaded db
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Get summary statistics for %s' % downloaded_pfam_db))
    hmmstat_cmd_db = 'hmmstat %s > %s' % (downloaded_pfam_db, downloaded_pfam_db_stats)
    os.system(hmmstat_cmd_db)

    # parse Pfam-A.hmm.stat and get updated id file
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Parse summary statistics for %s' % downloaded_pfam_db))
    id_to_version_dict_new = {}
    Pfam_41_ids_updated_handle = open(Pfam_41_ids_updated, 'w')
    for each_hmm in open(downloaded_pfam_db_stats):

        if not each_hmm.startswith('#') and (each_hmm.strip() != ''):

            each_hmm_split = each_hmm.strip().split(' ')

            # remove '' from each_hmm_split
            each_hmm_split_no_space = []
            for i in each_hmm_split:
                if i != '':
                    each_hmm_split_no_space.append(i)

            hmm_id_version = each_hmm_split_no_space[2]
            hmm_id_no_version = hmm_id_version.split('.')[0]
            id_to_version_dict_new[hmm_id_no_version] = hmm_id_version

            if hmm_id_no_version in Pfam_41_id_list_no_version:
                Pfam_41_ids_updated_handle.write('%s\n' % hmm_id_version)

    Pfam_41_ids_updated_handle.close()

    # extract updated hmm profile from downloaded db
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Extract newest Pfam profiles from %s' % downloaded_pfam_db))
    hmmfetch_cmd = 'hmmfetch -f %s %s > %s' % (downloaded_pfam_db, Pfam_41_ids_updated, Pfam_profiles_updated)
    os.system(hmmfetch_cmd)

    # convert TIGRFAM hmm profiles to HMMER3 ASCII format
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Convert needed TIGRFAM profiles to HMMER3 ASCII format'))
    hmmconvert_cmd_TIGR00344 = 'hmmconvert %s/TIGR00344.HMM > TIGR00344.HMM' % TIGRFAMs_profiles_folder
    hmmconvert_cmd_TIGR00422 = 'hmmconvert %s/TIGR00422.HMM > TIGR00422.HMM' % TIGRFAMs_profiles_folder
    os.system(hmmconvert_cmd_TIGR00344)
    os.system(hmmconvert_cmd_TIGR00422)

    # combine all 43 profiles together
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Combine Pfam and TIGRFAM profiles together'))
    cat_cmd = 'cat %s TIGR00344.HMM TIGR00422.HMM > %s' % (Pfam_profiles_updated, updated_43_profiles)
    os.system(cat_cmd)

    ######################################## get log file ########################################

    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Prepare log file'))
    update_log_file_handle = open(update_log_file, 'w')
    for each_id in Pfam_41_id_list_no_version:
        update_log_file_handle.write('%s --> %s\n' % (id_to_version_dict_old[each_id], id_to_version_dict_new[each_id]))
    update_log_file_handle.close()

    ######################################## delete tmp files ########################################

    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Delete temporary files'))
    os.system('rm %s' % current_43_hmm_file_stat)
    os.system('rm %s' % downloaded_pfam_db_stats)
    os.system('rm %s' % Pfam_41_ids_updated)
    os.system('rm %s' % Pfam_profiles_updated)
    os.system('rm TIGR00344.HMM')
    os.system('rm TIGR00422.HMM')

    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Done, updated profiles exported to %s' % Pfam_profiles_updated))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-hmm',  required=True,                help='MetaCHIP_phylo.hmm file')
    parser.add_argument('-p_db', required=False, default=None, help='Pfam db file, e.g. Pfam-A.hmm')
    parser.add_argument('-t_db', required=False, default=None, help='TIGRFAMs db folder, e.g. TIGRFAMs_14.0_HMM')

    args = vars(parser.parse_args())

    update_hmms(args)
