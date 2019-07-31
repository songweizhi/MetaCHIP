#!/usr/bin/env python
from __future__ import division
import os
import re
import glob
import shutil
import argparse
import warnings
from datetime import datetime
from Bio import SeqIO, AlignIO, Align
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio import SeqFeature as SF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import multiprocessing as mp


get_SCG_tree_usage = '''
===================================== get SCG tree example commands =====================================

# for completed genome
MetaCHIP get_SCG_tree -i genomes -p NorthSea -x fasta -t 4 -nonmeta

# for metagenome-assembled genomes (MAGs) 
MetaCHIP get_SCG_tree -i genomes -p NorthSea -x fasta -t 4

# Software dependencies:
Prodigal, HMMER, Mafft and FastTree

=========================================================================================================
'''


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


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


def remove_empty_element(list_in):

    list_out = []
    for each_element in list_in:
        if each_element != '':
            list_out.append(each_element)

    return list_out


def get_program_path_dict(pwd_cfg_file):
    program_path_dict = {}
    for each in open(pwd_cfg_file):
        each_split = each.strip().split('=')
        program_name = each_split[0]
        program_path = each_split[1]

        # remove space if there are
        if program_name[-1] == ' ':
            program_name = program_name[:-1]
        if program_path[0] == ' ':
            program_path = program_path[1:]
        program_path_dict[program_name] = program_path

    return program_path_dict


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def export_aa_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.protein)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def remove_low_cov_and_consensus_columns(alignment_file_in, minimal_cov, min_consensus, alignment_file_out):

    def list_to_segments(list_in):

        segments_out = []
        current_element = None
        current_segment = [None, None]
        for each_element in list_in:

            # for the first ellment
            if current_element == None:
                current_element = each_element
                current_segment = [each_element, each_element]

            elif each_element == current_element + 1:
                current_segment[1] = each_element
                current_element = each_element

            elif each_element != current_element + 1:

                # add segment to list
                segments_out.append(current_segment)

                # resetting segment
                current_segment = [each_element, each_element]
                current_element = each_element

        # add segment to list
        segments_out.append(current_segment)

        return segments_out

    def remove_columns_from_msa(alignment_in, cols_to_remove):

        # get 0 based index of all wanted columns
        cols_to_remove_0_base = [(i - 1) for i in cols_to_remove]
        aln_cols_index_all = list(range(alignment_in.get_alignment_length()))
        aln_cols_index_wanted = []
        for i in aln_cols_index_all:
            if i not in cols_to_remove_0_base:
                aln_cols_index_wanted.append(i)

        # get wanted alignment segments
        wanted_segments = list_to_segments(aln_cols_index_wanted)

        # create an empty Alignment object
        alignment_new = Align.MultipleSeqAlignment([])
        for sequence in alignment_in:
            new_seq_object = Seq('')
            new_seq_record = SeqRecord(new_seq_object)
            new_seq_record.id = sequence.id
            new_seq_record.description = sequence.description
            alignment_new.append(new_seq_record)

        # add wanted columns to empty Alignment object
        for segment in wanted_segments:

            # for single column segment
            if segment[0] == segment[1]:
                segment_value = alignment_in[:, segment[0]]

                m = 0
                for each_seq in alignment_new:
                    each_seq.seq = Seq(str(each_seq.seq) + segment_value[m])
                    m += 1

            # for multiple columns segment
            else:
                segment_value = alignment_in[:, (segment[0]):(segment[1] + 1)]
                alignment_new += segment_value

        return alignment_new

    def remove_low_cov_columns(alignment_in, min_cov_cutoff):

        # get columns with low coverage
        sequence_number = len(alignment_in)
        total_col_num = alignment_in.get_alignment_length()
        low_cov_columns = []
        n = 0
        while n < total_col_num:
            current_column = alignment_in[:, n]
            dash_number = current_column.count('-')
            gap_percent = (dash_number / sequence_number) * 100

            if gap_percent > min_cov_cutoff:
                low_cov_columns.append(n + 1)

            n += 1

        # remove identified columns
        alignment_new = remove_columns_from_msa(alignment_in, low_cov_columns)

        return alignment_new

    def remove_low_consensus_columns(alignment_in, min_css_cutoff):

        # get columns with low coverage
        sequence_number = len(alignment_in)
        total_col_num = alignment_in.get_alignment_length()
        low_css_columns = []
        n = 0
        while n < total_col_num:
            current_column = alignment_in[:, n]

            # get all aa in current column
            aa_list = set()
            for aa in current_column:
                aa_list.add(aa)

            # get maximum aa percent
            most_abundant_aa_percent = 0
            for each_aa in aa_list:
                each_aa_percent = (current_column.count(each_aa) / sequence_number) * 100
                if each_aa_percent > most_abundant_aa_percent:
                    most_abundant_aa_percent = each_aa_percent

            # if maximum percent lower than provided cutoff, add current column to low consensus column list
            if most_abundant_aa_percent < min_css_cutoff:
                low_css_columns.append(n + 1)

            n += 1

        # remove identified columns
        alignment_new = remove_columns_from_msa(alignment_in, low_css_columns)

        return alignment_new

    # read in alignment
    alignment = AlignIO.read(alignment_file_in, "fasta")

    # remove_low_cov_columns
    alignment_cov = remove_low_cov_columns(alignment, minimal_cov)

    # remove_low_consensus_columns
    alignment_cov_css = remove_low_consensus_columns(alignment_cov, min_consensus)

    # write filtered alignment
    alignment_file_out_handle = open(alignment_file_out, 'w')
    for each_seq in alignment_cov_css:
        alignment_file_out_handle.write('>%s\n' % str(each_seq.id))
        alignment_file_out_handle.write('%s\n' % str(each_seq.seq))
    alignment_file_out_handle.close()


def prodigal_parser(seq_file, sco_file, prefix, output_folder):

    bin_ffn_file =     '%s.ffn' % prefix
    bin_faa_file =     '%s.faa' % prefix
    pwd_bin_ffn_file = '%s/%s'  % (output_folder, bin_ffn_file)
    pwd_bin_faa_file = '%s/%s'  % (output_folder, bin_faa_file)

    # get sequence id list
    id_to_sequence_dict = {}
    sequence_id_list = []
    for each_seq in SeqIO.parse(seq_file, 'fasta'):
        id_to_sequence_dict[each_seq.id] = str(each_seq.seq)
        sequence_id_list.append(each_seq.id)


    # get sequence to cds dict and sequence to transl_table dict
    current_seq_id = ''
    current_transl_table = ''
    current_seq_csd_list = []
    seq_to_cds_dict = {}
    seq_to_transl_table_dict = {}
    for each_cds in open(sco_file):
        if each_cds.startswith('# Sequence Data'):

            # add to dict
            if current_seq_id != '':
                seq_to_cds_dict[current_seq_id] = current_seq_csd_list
                seq_to_transl_table_dict[current_seq_id] = current_transl_table

            # reset value
            current_seq_id = each_cds.strip().split('=')[-1][1:-1].split(' ')[0]
            current_transl_table = ''
            current_seq_csd_list = []

        elif each_cds.startswith('# Model Data'):
            current_transl_table = each_cds.strip().split(';')[-2].split('=')[-1]

        else:
            current_seq_csd_list.append('_'.join(each_cds.strip().split('_')[1:]))

    seq_to_cds_dict[current_seq_id] = current_seq_csd_list
    seq_to_transl_table_dict[current_seq_id] = current_transl_table


    bin_ffn_file_handle = open(pwd_bin_ffn_file, 'w')
    bin_faa_file_handle = open(pwd_bin_faa_file, 'w')
    gene_index = 1
    for seq_id in sequence_id_list:

        # create SeqRecord
        current_sequence = Seq(id_to_sequence_dict[seq_id])
        current_SeqRecord = SeqRecord(current_sequence, id=seq_id)
        current_SeqRecord.seq.alphabet = generic_dna
        transl_table = seq_to_transl_table_dict[seq_id]

        # add SeqFeature to SeqRecord
        for cds in seq_to_cds_dict[seq_id]:

            # define locus_tag id
            locus_tag_id = '%s_%s' % (prefix, "{:0>5}".format(gene_index))

            # define FeatureLocation
            cds_split = cds.split('_')
            cds_start = SF.ExactPosition(int(cds_split[0]))
            cds_end = SF.ExactPosition(int(cds_split[1]))
            cds_strand = cds_split[2]
            current_strand = None
            if cds_strand == '+':
                current_strand = 1
            if cds_strand == '-':
                current_strand = -1
            current_feature_location = FeatureLocation(cds_start, cds_end, strand=current_strand)

            # get nc sequence
            sequence_nc = ''
            if cds_strand == '+':
                sequence_nc = id_to_sequence_dict[seq_id][cds_start-1:cds_end]
            if cds_strand == '-':
                sequence_nc = str(Seq(id_to_sequence_dict[seq_id][cds_start-1:cds_end], generic_dna).reverse_complement())

            # translate to aa sequence
            sequence_aa = str(SeqRecord(Seq(sequence_nc)).seq.translate(table=transl_table))

            # remove * at the end
            sequence_aa = sequence_aa[:-1]

            # export nc and aa sequences
            export_dna_record(sequence_nc, locus_tag_id, '', bin_ffn_file_handle)
            export_aa_record(sequence_aa, locus_tag_id, '', bin_faa_file_handle)

            # Define feature type
            current_feature_type = 'CDS'

            # Define feature qualifiers
            current_qualifiers_dict = {}
            current_qualifiers_dict['locus_tag'] = locus_tag_id
            current_qualifiers_dict['transl_table'] = transl_table
            current_qualifiers_dict['translation'] = sequence_aa

            # Create a SeqFeature
            current_feature = SeqFeature(current_feature_location, type=current_feature_type, qualifiers=current_qualifiers_dict)

            # Append Feature to SeqRecord
            current_SeqRecord.features.append(current_feature)
            gene_index += 1

    bin_ffn_file_handle.close()
    bin_faa_file_handle.close()


def sep_combined_hmm(combined_hmm_file, hmm_profile_sep_folder, hmmfetch_exe, pwd_hmmstat_exe):

    # extract hmm profile id from phylo.hmm
    pwd_phylo_hmm_stat_txt = '%s/phylo.hmm.stat.txt' % hmm_profile_sep_folder
    hmmstat_cmd = '%s %s > %s' % (pwd_hmmstat_exe, combined_hmm_file, pwd_phylo_hmm_stat_txt)
    os.system(hmmstat_cmd)

    # get hmm profile id file
    hmm_id_list = []
    for each_profile in open(pwd_phylo_hmm_stat_txt):
        if not each_profile.startswith('#'):
            each_profile_split = each_profile.strip().split(' ')
            if each_profile_split != ['']:
                each_profile_split_no_space = []
                for each_element in each_profile_split:
                    if each_element != '':
                        each_profile_split_no_space.append(each_element)
                hmm_id_list.append(each_profile_split_no_space[2])

    for each_hmm_id in hmm_id_list:
        hmmfetch_cmd = '%s %s %s > %s/%s.hmm' % (hmmfetch_exe, combined_hmm_file, each_hmm_id, hmm_profile_sep_folder, each_hmm_id)
        os.system(hmmfetch_cmd)


def prodigal_worker(argument_list):

    input_genome = argument_list[0]
    input_genome_folder = argument_list[1]
    pwd_prodigal_exe = argument_list[2]
    nonmeta_mode = argument_list[3]
    pwd_prodigal_output_folder = argument_list[4]

    # prepare command (according to Prokka)
    input_genome_basename, input_genome_ext = os.path.splitext(input_genome)
    pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    pwd_output_sco = '%s/%s.sco' % (pwd_prodigal_output_folder, input_genome_basename)

    prodigal_cmd_meta = '%s -f sco -q -c -m -g 11 -p meta -i %s -o %s' % (
    pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)
    prodigal_cmd_nonmeta = '%s -f sco -q -c -m -g 11 -i %s -o %s' % (
    pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)

    if nonmeta_mode is True:
        prodigal_cmd = prodigal_cmd_nonmeta
    else:
        prodigal_cmd = prodigal_cmd_meta

    os.system(prodigal_cmd)

    # prepare ffn, faa and gbk files from prodigal output
    prodigal_parser(pwd_input_genome, pwd_output_sco, input_genome_basename, pwd_prodigal_output_folder)


def hmmsearch_worker(argument_list):

    faa_file_basename = argument_list[0]
    pwd_SCG_tree_wd = argument_list[1]
    pwd_hmmsearch_exe = argument_list[2]
    path_to_hmm = argument_list[3]
    pwd_faa_folder = argument_list[4]

    # run hmmsearch
    pwd_faa_file = '%s/%s.faa' % (pwd_faa_folder, faa_file_basename)
    os.system('%s -o /dev/null --domtblout %s/%s_hmmout.tbl %s %s' % (pwd_hmmsearch_exe, pwd_SCG_tree_wd, faa_file_basename, path_to_hmm, pwd_faa_file))

    # Reading the protein file in a dictionary
    proteinSequence = {}
    for seq_record in SeqIO.parse(pwd_faa_file, 'fasta'):
        proteinSequence[seq_record.id] = str(seq_record.seq)

    # Reading the hmmersearch table/extracting the protein part found beu hmmsearch out of the protein/Writing
    # each protein sequence that was extracted to a fasta file (one for each hmm in phylo.hmm
    hmm_id = ''
    hmm_name = ''
    hmm_pos1 = 0
    hmm_pos2 = 0
    hmm_score = 0
    pwd_hmmout_tbl = pwd_SCG_tree_wd + '/' + faa_file_basename + '_hmmout.tbl'
    with open(pwd_hmmout_tbl, 'r') as tbl:
        for line in tbl:
            if line[0] == "#": continue
            line = re.sub('\s+', ' ', line)
            splitLine = line.split(' ')

            if (hmm_id == ''):
                hmm_id = splitLine[4]
                hmm_name = splitLine[0]
                hmm_pos1 = int(splitLine[17]) - 1
                hmm_pos2 = int(splitLine[18])
                hmm_score = float(splitLine[13])
            elif (hmm_id == splitLine[4]):
                if (float(splitLine[13]) > hmm_score):
                    hmm_name = splitLine[0]
                    hmm_pos1 = int(splitLine[17]) - 1
                    hmm_pos2 = int(splitLine[18])
                    hmm_score = float(splitLine[13])
            else:
                file_out = open(pwd_SCG_tree_wd + '/' + hmm_id + '.fasta', 'a+')
                file_out.write('>' + faa_file_basename + '\n')
                if hmm_name != '':
                    seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
                    file_out.write(str(seq) + '\n')
                    file_out.close()
                hmm_id = splitLine[4]
                hmm_name = splitLine[0]
                hmm_pos1 = int(splitLine[17]) - 1
                hmm_pos2 = int(splitLine[18])
                hmm_score = float(splitLine[13])

        else:
            file_out = open(pwd_SCG_tree_wd + '/' + hmm_id + '.fasta', 'a+')
            file_out.write('>' + faa_file_basename + '\n')
            if hmm_name != '':
                seq = str(proteinSequence[hmm_name][hmm_pos1:hmm_pos2])
                file_out.write(str(seq) + '\n')
                file_out.close()


def convert_hmmalign_output(align_in, align_out):

    # read in alignment
    sequence_id_list = []
    sequence_seq_dict = {}
    for aligned_seq in open(align_in):
        aligned_seq_split = aligned_seq.strip().split(' ')
        aligned_seq_split = remove_empty_element(aligned_seq_split)

        if aligned_seq_split != []:
            aligned_seq_id = aligned_seq_split[0]
            aligned_seq_seq = aligned_seq_split[1]

            # add id to sequence id list
            if aligned_seq_id not in sequence_id_list:
                sequence_id_list.append(aligned_seq_id)

            # add seq to sequence seq dict
            if aligned_seq_id not in sequence_seq_dict:
                sequence_seq_dict[aligned_seq_id] = aligned_seq_seq
            else:
                sequence_seq_dict[aligned_seq_id] += aligned_seq_seq

    # write out
    align_out_handle = open(align_out, 'w')
    for sequence_id in sequence_id_list:
        sequence_seq = sequence_seq_dict[sequence_id]
        align_out_handle.write('>%s\n' % sequence_id)
        align_out_handle.write('%s\n' % sequence_seq)
    align_out_handle.close()


def hmmalign_worker(argument_list):
    fastaFile_basename = argument_list[0]
    pwd_SCG_tree_wd = argument_list[1]
    pwd_hmm_profile_folder = argument_list[2]
    pwd_hmmalign_exe = argument_list[3]

    pwd_hmm_file =    '%s/%s.hmm'               % (pwd_hmm_profile_folder, fastaFile_basename)
    pwd_seq_in =      '%s/%s.fasta'             % (pwd_SCG_tree_wd, fastaFile_basename)
    pwd_aln_out_tmp = '%s/%s_aligned_tmp.fasta' % (pwd_SCG_tree_wd, fastaFile_basename)
    pwd_aln_out =     '%s/%s_aligned.fasta'     % (pwd_SCG_tree_wd, fastaFile_basename)

    hmmalign_cmd = '%s --trim --outformat PSIBLAST %s %s > %s ; rm %s' % (pwd_hmmalign_exe, pwd_hmm_file, pwd_seq_in, pwd_aln_out_tmp, pwd_seq_in)
    os.system(hmmalign_cmd)

    # convert alignment format
    convert_hmmalign_output(pwd_aln_out_tmp, pwd_aln_out)

    # remove tmp alignment
    os.system('rm %s' % pwd_aln_out_tmp)


config_file_path = '/Users/songweizhi/PycharmProjects/MyBioTools/MyBioTools'

config_dict = {'prodigal'         : 'prodigal',
               'hmmsearch'        : 'hmmsearch',
               'hmmfetch'         : 'hmmfetch',
               'hmmalign'         : 'hmmalign',
               'hmmstat'          : 'hmmstat',
               'mafft'            : 'mafft',
               'blastp'           : 'blastp',
               'blastn'           : 'blastn',
               'makeblastdb'      : 'makeblastdb',
               'fasttree'         : 'FastTree',
               'ranger_mac'       : '%s/Ranger-DTL-Dated.mac'   % config_file_path,  # do not edit this line
               'ranger_linux'     : '%s/Ranger-DTL-Dated.linux' % config_file_path,  # do not edit this line
               'path_to_hmm'      : '%s/MetaCHIP_phylo.hmm'     % config_file_path,  # do not edit this line
               'circos_HGT_R'     : '%s/MetaCHIP_circos_HGT.R'  % config_file_path,   # do not edit this line
               'cdd2cog_perl'     : '%s/cdd2cog.pl' % config_file_path,
               'get_sankey_plot_R': '%s/get_sankey_plot.R' % config_file_path
               }


def get_SCG_tree(args, config_dict):

    # read in arguments
    input_genome_folder =   args['i']
    output_prefix =         args['p']
    file_extension =        args['x']
    num_threads =           args['t']
    nonmeta_mode =          args['nonmeta']

    # read in config file
    path_to_hmm =           config_dict['path_to_hmm']
    pwd_prodigal_exe =      config_dict['prodigal']
    pwd_hmmsearch_exe =     config_dict['hmmsearch']
    pwd_hmmfetch_exe =      config_dict['hmmfetch']
    pwd_hmmalign_exe =      config_dict['hmmalign']
    pwd_hmmstat_exe =       config_dict['hmmstat']
    pwd_fasttree_exe =      config_dict['fasttree']

    warnings.filterwarnings("ignore")
    minimal_cov_in_msa = 50
    min_consensus_in_msa = 25
    keep_quiet = False


    #################################################### check input ###################################################

    # check whether input genome exist
    input_genome_file_re = '%s/*.%s' % (input_genome_folder, file_extension)
    input_genome_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_file_re)]
    if input_genome_file_name_list == []:
        print('No input genome detected, program exited!')
        exit()


    ############################################# define file/folder names #############################################

    get_SCG_tree_wd =                    '%s_get_SCG_tree_wd'                   % (output_prefix)
    prodigal_output_folder =             '%s_1_prodigal_output'                 % (output_prefix)
    extract_and_align_SCG_wd =           '%s_2_extract_and_align_SCGs'          % (output_prefix)
    combined_alignment_file_tmp =        '%s_SCG_tree.aln'                      % (output_prefix)
    combined_alignment_file =            '%s_SCG_tree_cov%s_css%s.aln'          % (output_prefix, minimal_cov_in_msa, min_consensus_in_msa)
    newick_tree_file =                   '%s_SCG_tree.newick'                   % (output_prefix)
    hmm_profile_sep_folder =             '%s_hmm_profile_fetched'               % (output_prefix)

    pwd_log_file =                       '%s/%s_get_SCG_tree.log'               % (get_SCG_tree_wd, output_prefix)
    pwd_prodigal_output_folder =         '%s/%s'                                % (get_SCG_tree_wd, prodigal_output_folder)
    pwd_extract_and_align_SCG_wd =       '%s/%s'                                % (get_SCG_tree_wd, extract_and_align_SCG_wd)
    pwd_combined_alignment_file_tmp =    '%s/%s'                                % (get_SCG_tree_wd, combined_alignment_file_tmp)
    pwd_combined_alignment_file =        '%s/%s'                                % (get_SCG_tree_wd, combined_alignment_file)
    pwd_hmm_profile_sep_folder =         '%s/%s/%s'                             % (get_SCG_tree_wd, extract_and_align_SCG_wd, hmm_profile_sep_folder)
    pwd_newick_tree_file =               '%s/%s'                                % (get_SCG_tree_wd, newick_tree_file)


    # create wd
    force_create_folder(get_SCG_tree_wd)


    ######################################## run prodigal with multiprocessing #########################################

    # for report and log
    report_and_log(('Running Prodigal with %s cores for input genomes' % num_threads), pwd_log_file, keep_quiet)

    # create prodigal output folder
    force_create_folder(pwd_prodigal_output_folder)

    # get input genome list
    input_genome_file_re = '%s/*.%s' % (input_genome_folder, file_extension)
    input_genome_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_file_re)]

    # prepare arguments for prodigal_worker
    list_for_multiple_arguments_Prodigal = []
    for input_genome in input_genome_file_name_list:
        list_for_multiple_arguments_Prodigal.append([input_genome, input_genome_folder, pwd_prodigal_exe, nonmeta_mode, pwd_prodigal_output_folder])

    # run prodigal with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(prodigal_worker, list_for_multiple_arguments_Prodigal)
    pool.close()
    pool.join()


    ########################################### get species tree (hmmsearch) ###########################################

    # create wd
    force_create_folder(pwd_extract_and_align_SCG_wd)

    # for report and log
    report_and_log(('Running Hmmsearch with %s cores' % num_threads), pwd_log_file, keep_quiet)

    faa_file_re = '%s/*.faa' % pwd_prodigal_output_folder
    faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
    faa_file_list = sorted(faa_file_list)

    faa_file_basename_list = []
    for faa_file in faa_file_list:
        faa_file_basename, faa_file_extension = os.path.splitext(faa_file)
        faa_file_basename_list.append(faa_file_basename)

    # prepare arguments for hmmsearch_worker
    list_for_multiple_arguments_hmmsearch = []
    for faa_file_basename in faa_file_basename_list:
        list_for_multiple_arguments_hmmsearch.append([faa_file_basename, pwd_extract_and_align_SCG_wd, pwd_hmmsearch_exe, path_to_hmm, pwd_prodigal_output_folder])

    # run hmmsearch with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(hmmsearch_worker, list_for_multiple_arguments_hmmsearch)
    pool.close()
    pool.join()


    ############################################# get species tree (hmmalign) #############################################

    # for report and log
    report_and_log(('Running Hmmalign with %s cores' % num_threads), pwd_log_file, keep_quiet)

    # fetch combined hmm profiles
    force_create_folder(pwd_hmm_profile_sep_folder)
    sep_combined_hmm(path_to_hmm, pwd_hmm_profile_sep_folder, pwd_hmmfetch_exe, pwd_hmmstat_exe)

    # Call hmmalign to align all single fasta files with hmms
    files = os.listdir(pwd_extract_and_align_SCG_wd)
    fastaFiles = [i for i in files if i.endswith('.fasta')]

    # prepare arguments for hmmalign_worker
    list_for_multiple_arguments_hmmalign = []
    for fastaFile in fastaFiles:

        fastaFiles_basename = '.'.join(fastaFile.split('.')[:-1])
        list_for_multiple_arguments_hmmalign.append([fastaFiles_basename, pwd_extract_and_align_SCG_wd, pwd_hmm_profile_sep_folder, pwd_hmmalign_exe])

    # run hmmalign with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(hmmalign_worker, list_for_multiple_arguments_hmmalign)
    pool.close()
    pool.join()


    ################################### get species tree (Concatenating alignments) ####################################

    # for report and log
    report_and_log('Concatenating alignments', pwd_log_file, keep_quiet)

    # concatenating the single alignments
    concatAlignment = {}
    for element in faa_file_basename_list:
        concatAlignment[element] = ''

    # Reading all single alignment files and append them to the concatenated alignment
    files = os.listdir(pwd_extract_and_align_SCG_wd)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    for faa_file_basename in fastaFiles:
        fastaFile = pwd_extract_and_align_SCG_wd + '/' + faa_file_basename
        proteinSequence = {}
        alignmentLength = 0
        for seq_record_2 in SeqIO.parse(fastaFile, 'fasta'):
            proteinName = seq_record_2.id
            proteinSequence[proteinName] = str(seq_record_2.seq)
            alignmentLength = len(proteinSequence[proteinName])

        for element in faa_file_basename_list:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    file_out = open(pwd_combined_alignment_file_tmp, 'w')
    for element in faa_file_basename_list:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()

    # remove columns with low coverage and low consensus
    report_and_log(('Removing columns from concatenated alignment represented by <%s%s of genomes and with an amino acid consensus <%s%s' % (minimal_cov_in_msa, '%', min_consensus_in_msa, '%')), pwd_log_file, keep_quiet)
    remove_low_cov_and_consensus_columns(pwd_combined_alignment_file_tmp, minimal_cov_in_msa, min_consensus_in_msa, pwd_combined_alignment_file)


    ########################################### get species tree (fasttree) ############################################

    # for report and log
    report_and_log('Running FastTree', pwd_log_file, keep_quiet)

    # calling fasttree for tree calculation
    fasttree_cmd = '%s -quiet %s > %s' % (pwd_fasttree_exe, pwd_combined_alignment_file, pwd_newick_tree_file)
    os.system(fasttree_cmd)

    # for report and log
    report_and_log(('SCG tree exported to: %s' % newick_tree_file), pwd_log_file, keep_quiet)


    ############################################## remove temporary files ##############################################

    # remove temporary files
    report_and_log(('Deleting temporary files'), pwd_log_file, keep_quiet)

    os.system('rm -r %s' % pwd_combined_alignment_file_tmp)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    # arguments for PI
    parser.add_argument('-i',             required=True,  help='input genome folder')
    parser.add_argument('-p',             required=True,  help='output prefix')
    parser.add_argument('-x',             required=False, default='fasta',     help='file extension')
    parser.add_argument('-nonmeta',       required=False, action="store_true", help='annotate Non-metagenome-assembled genomes (Non-MAGs)')
    parser.add_argument('-t',             required=False, type=int, default=1, help='number of threads, default: 1')

    args = vars(parser.parse_args())

    get_SCG_tree(args, config_dict)

