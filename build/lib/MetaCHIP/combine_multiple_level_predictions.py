import os
import copy
import glob
import shutil
import numpy as np
import pandas as pd
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from MetaCHIP.MetaCHIP_config import config_dict


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


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def combine_PG_output(PG_output_file_list_with_path, output_prefix, detection_ranks, combined_PG_output_normal):

    HGT_identity_dict = {}
    HGT_end_match_dict = {}
    HGT_full_length_match_dict = {}
    HGT_direction_dict = {}
    HGT_occurence_dict = {}
    HGT_concatenated_list = []
    for pwd_PG_output_file in PG_output_file_list_with_path:
        file_path, file_name = os.path.split(pwd_PG_output_file)
        taxon_rank = file_name[len(output_prefix) + 1]
        if taxon_rank in detection_ranks:
            for PG_HGT in open(pwd_PG_output_file):
                if not PG_HGT.startswith('Gene_1'):
                    PG_HGT_split = PG_HGT.strip().split('\t')

                    gene_1 = PG_HGT_split[0]
                    gene_2 = PG_HGT_split[1]
                    identity = float(PG_HGT_split[4])
                    end_match = PG_HGT_split[5]
                    full_length_match = PG_HGT_split[6]
                    direction = PG_HGT_split[7]
                    concatenated = '%s___%s' % (gene_1, gene_2)

                    if concatenated not in HGT_concatenated_list:
                        HGT_concatenated_list.append(concatenated)

                    # store in dict
                    if concatenated not in HGT_identity_dict:
                        HGT_identity_dict[concatenated] = identity

                    if concatenated not in HGT_end_match_dict:
                        HGT_end_match_dict[concatenated] = end_match

                    if concatenated not in HGT_full_length_match_dict:
                        HGT_full_length_match_dict[concatenated] = full_length_match

                    if direction != 'NA':
                        if concatenated not in HGT_direction_dict:
                            HGT_direction_dict[concatenated] = [direction]
                        else:
                            HGT_direction_dict[concatenated].append(direction)

                    if direction != 'NA':
                        if concatenated not in HGT_occurence_dict:
                            HGT_occurence_dict[concatenated] = [taxon_rank]
                        else:
                            HGT_occurence_dict[concatenated].append(taxon_rank)


    detection_ranks_all = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    detection_ranks_list = []
    for each_rank in detection_ranks_all:
        if each_rank in detection_ranks:
            detection_ranks_list.append(each_rank)


    HGT_occurence_dict_0_1_format = {}
    for each_HGT in HGT_occurence_dict:
        occurence_str = ''
        for each_level in detection_ranks_list:
            if each_level in HGT_occurence_dict[each_HGT]:
                occurence_str += '1'
            else:
                occurence_str += '0'
        HGT_occurence_dict_0_1_format[each_HGT] = occurence_str


    combined_output_handle_normal = open(combined_PG_output_normal, 'w')

    combined_output_handle_normal.write('Gene_1\tGene_2\tIdentity\toccurence(%s)\tend_match\tfull_length_match\tdirection\n' % detection_ranks)

    for concatenated_HGT in sorted(HGT_concatenated_list):
        concatenated_HGT_split = concatenated_HGT.split('___')

        concatenated_HGT_direction = 'NA'
        if concatenated_HGT in HGT_direction_dict:

            concatenated_HGT_direction_list = HGT_direction_dict[concatenated_HGT]
            concatenated_HGT_direction_list_uniq = unique_list_elements(concatenated_HGT_direction_list)

            if len(concatenated_HGT_direction_list_uniq) == 1:
                concatenated_HGT_direction = concatenated_HGT_direction_list[0]
            else:
                concatenated_HGT_direction = 'both'
                for HGT_direction in concatenated_HGT_direction_list_uniq:
                    HGT_direction_freq = (concatenated_HGT_direction_list.count(HGT_direction))/len(concatenated_HGT_direction_list)
                    if HGT_direction_freq > 0.5:
                        concatenated_HGT_direction = HGT_direction + '(' + str(float("{0:.2f}".format(HGT_direction_freq*100))) + '%)'

        if concatenated_HGT in HGT_occurence_dict_0_1_format:
            occurence_formatted = HGT_occurence_dict_0_1_format[concatenated_HGT]
        else:
            occurence_formatted = '0'*len(detection_ranks_list)

        for_out = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (concatenated_HGT_split[0],
                                                    concatenated_HGT_split[1],
                                                    HGT_identity_dict[concatenated_HGT],
                                                    occurence_formatted,
                                                    HGT_end_match_dict[concatenated_HGT],
                                                    HGT_full_length_match_dict[concatenated_HGT],
                                                    concatenated_HGT_direction)

        if (HGT_end_match_dict[concatenated_HGT] == 'no') and (HGT_full_length_match_dict[concatenated_HGT] == 'no') and (concatenated_HGT_direction != 'NA') and (concatenated_HGT_direction != 'both'):
            combined_output_handle_normal.write(for_out)

    combined_output_handle_normal.close()


def combine_multiple_level_predictions(args, config_dict):

    output_prefix =             args['p']
    grouping_level =            args['r']
    grouping_file =             args['g']
    cover_cutoff =              args['cov']
    align_len_cutoff =          args['al']
    flanking_length_kbp =       args['flk']
    identity_percentile =       args['ip']
    end_match_identity_cutoff = args['ei']
    keep_quiet =                args['quiet']
    keep_temp =                 args['tmp']

    detection_rank_list = args['r']

    # for single level detection
    if len(detection_rank_list) == 1:

        pwd_MetaCHIP_op_folder_re = '%s_MetaCHIP_wd/%s_%s*_HGTs_ip%s_al%sbp_c%s_ei%sbp_f%skbp' % (output_prefix, output_prefix, grouping_level, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)
        pwd_MetaCHIP_op_folder = [os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)][0]

        print(pwd_MetaCHIP_op_folder)


    # for multiple level detection
    if len(detection_rank_list) > 1:

        pwd_detected_HGT_txt_list = []
        for detection_rank in detection_rank_list:

            pwd_MetaCHIP_op_folder_re = '%s_MetaCHIP_wd/%s_%s*_HGTs_ip%s_al%sbp_c%s_ei%sbp_f%skbp' % (output_prefix, output_prefix, detection_rank, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)
            MetaCHIP_op_folder = [os.path.basename(file_name) for file_name in glob.glob(pwd_MetaCHIP_op_folder_re)][0]

            group_num = int(MetaCHIP_op_folder[len(output_prefix)+1:].split('_')[0][1:])
            pwd_detected_HGT_txt = '%s_MetaCHIP_wd/%s/%s_%s%s_HGTs_PG.txt' % (output_prefix, MetaCHIP_op_folder, output_prefix, detection_rank, group_num)

            pwd_detected_HGT_txt_list.append(pwd_detected_HGT_txt)

        print(pwd_detected_HGT_txt_list)

        pwd_combined_prediction_folder = '%s_MetaCHIP_wd/%s_combined_%s_HGTs_ip%s_al%sbp_c%s_ei%sbp_f%skbp' % (output_prefix, output_prefix, detection_rank_list, str(identity_percentile), str(align_len_cutoff), str(cover_cutoff), str(end_match_identity_cutoff), flanking_length_kbp)

        force_create_folder(pwd_combined_prediction_folder)

        genome_size_file =                          '%s_MetaCHIP_wd/%s_all_genome_size.txt'         % (output_prefix, output_prefix)
        pwd_detected_HGT_txt_combined =             '%s/%s_%s_detected_HGTs.txt'                    % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
        pwd_recipient_gene_seq_ffn =                '%s/%s_%s_detected_HGTs_recipient_genes.ffn'    % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
        pwd_recipient_gene_seq_faa =                '%s/%s_%s_detected_HGTs_recipient_genes.faa'    % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
        pwd_donor_gene_seq_ffn =                    '%s/%s_%s_detected_HGTs_donor_genes.ffn'        % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)
        pwd_donor_gene_seq_faa =                    '%s/%s_%s_detected_HGTs_donor_genes.faa'        % (pwd_combined_prediction_folder, output_prefix, detection_rank_list)

        # combine prediction
        combine_PG_output(pwd_detected_HGT_txt_list, output_prefix, detection_rank_list, pwd_detected_HGT_txt_combined)


        ############################################### extract sequences ##############################################

        # get recipient and donor gene list
        recipient_gene_list = set()
        recipient_genome_list = []
        donor_gene_list = set()
        for each in open(pwd_detected_HGT_txt_combined):
            if not each.startswith('Gene_1'):

                each_split = each.strip().split('\t')
                gene_1 = each_split[0]
                gene_2 = each_split[1]
                gene_1_genome = '_'.join(gene_1.split('_')[:-1])
                gene_2_genome = '_'.join(gene_2.split('_')[:-1])
                direction = each_split[6]

                recipient_genome = direction.split('-->')[1]
                if '%)' in recipient_genome:
                    recipient_genome = recipient_genome.split('(')[0]
                recipient_genome_list.append(recipient_genome)

                if gene_1_genome == recipient_genome:
                    recipient_gene_list.add(gene_1)
                    donor_gene_list.add(gene_2)
                if gene_2_genome == recipient_genome:
                    recipient_gene_list.add(gene_2)
                    donor_gene_list.add(gene_1)


        # cat ffn and faa files from prodigal output folder
        pwd_combined_ffn = '%s/combined.ffn' % pwd_combined_prediction_folder
        os.system('cat %s_MetaCHIP_wd/%s_all_prodigal_output/*.ffn > %s' % (output_prefix, output_prefix, pwd_combined_ffn))

        pwd_recipient_gene_seq_ffn_handle = open(pwd_recipient_gene_seq_ffn, 'w')
        pwd_recipient_gene_seq_faa_handle = open(pwd_recipient_gene_seq_faa, 'w')
        pwd_donor_gene_seq_ffn_handle = open(pwd_donor_gene_seq_ffn, 'w')
        pwd_donor_gene_seq_faa_handle = open(pwd_donor_gene_seq_faa, 'w')

        for each_seq in SeqIO.parse(pwd_combined_ffn, 'fasta'):

            if str(each_seq.id) in recipient_gene_list:

                # write out nc sequences
                SeqIO.write(each_seq, pwd_recipient_gene_seq_ffn_handle, 'fasta')

                # write out aa sequences
                each_seq_aa = each_seq
                each_seq_aa.seq = each_seq_aa.seq.translate()
                SeqIO.write(each_seq_aa, pwd_recipient_gene_seq_faa_handle, 'fasta')

            if str(each_seq.id) in donor_gene_list:

                # write out nc sequences
                SeqIO.write(each_seq, pwd_donor_gene_seq_ffn_handle, 'fasta')

                # write out aa sequences
                each_seq_aa = each_seq
                each_seq_aa.seq = each_seq_aa.seq.translate()
                SeqIO.write(each_seq_aa, pwd_donor_gene_seq_faa_handle, 'fasta')


        pwd_recipient_gene_seq_ffn_handle.close()
        pwd_recipient_gene_seq_faa_handle.close()
        pwd_donor_gene_seq_ffn_handle.close()
        pwd_donor_gene_seq_faa_handle.close()

        os.remove(pwd_combined_ffn)


        ###################################### plot the number of HGT per genome #######################################


        #



os.chdir('/Users/songweizhi/Desktop/MetaCHIP_test')
args = {'p': 'NorthSea', 'r': 'pcof', 'g': None, 'cov': 75, 'al': 200, 'flk': 10, 'ip': 90, 'ei': 90, 't': 4, 'plot_iden': False, 'NoEbCheck': False, 'force': False, 'quiet': False, 'tmp': False}




combine_multiple_level_predictions(args, config_dict)





