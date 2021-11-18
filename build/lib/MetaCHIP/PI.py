#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com or t.thomas@unsw.edu.au

# MetaCHIP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MetaCHIP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import re
import glob
import shutil
import argparse
import warnings
import itertools
from time import sleep
from datetime import datetime
from string import ascii_uppercase
from Bio import SeqIO, AlignIO, Align
from Bio.Seq import Seq
from Bio import SeqFeature as SF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import multiprocessing as mp
from MetaCHIP.MetaCHIP_config import config_dict
from distutils.spawn import find_executable
warnings.filterwarnings("ignore")


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


def force_create_folder(folder_to_create):

    rm_rd = 0
    while os.path.isdir(folder_to_create) is True:
        shutil.rmtree(folder_to_create, ignore_errors=True)

        if rm_rd >= 10:
            print('Failed in removing %s, program exited!' % folder_to_create)
            exit()

        rm_rd += 1
        sleep(1)

    os.mkdir(folder_to_create)


def get_group_index_list():

    def iter_all_strings():
        size = 1
        while True:
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)
            size += 1

    group_index_list = []
    for s in iter_all_strings():
        group_index_list.append(s)
        if s == 'ZZZ':
            break

    return group_index_list


def prodigal_parser(seq_file, sco_file, prefix, output_folder):

    bin_ffn_file =     '%s.ffn' % prefix
    bin_faa_file =     '%s.faa' % prefix
    bin_gbk_file =     '%s.gbk' % prefix
    pwd_bin_ffn_file = '%s/%s'  % (output_folder, bin_ffn_file)
    pwd_bin_faa_file = '%s/%s'  % (output_folder, bin_faa_file)
    pwd_bin_gbk_file = '%s/%s'  % (output_folder, bin_gbk_file)

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
            current_seq_id = each_cds.strip().split(';seqhdr=')[1][1:-1].split(' ')[0]
            current_transl_table = ''
            current_seq_csd_list = []

        elif each_cds.startswith('# Model Data'):
            current_transl_table = each_cds.strip().split(';')[-2].split('=')[-1]

        else:
            current_seq_csd_list.append('_'.join(each_cds.strip().split('_')[1:]))

    seq_to_cds_dict[current_seq_id] = current_seq_csd_list
    seq_to_transl_table_dict[current_seq_id] = current_transl_table


    bin_gbk_file_handle = open(pwd_bin_gbk_file, 'w')
    bin_ffn_file_handle = open(pwd_bin_ffn_file, 'w')
    bin_faa_file_handle = open(pwd_bin_faa_file, 'w')
    gene_index = 1
    for seq_id in sequence_id_list:

        # create SeqRecord
        current_SeqRecord = SeqRecord(Seq(id_to_sequence_dict[seq_id]), id=seq_id, annotations={"molecule_type": "DNA"})
        transl_table = seq_to_transl_table_dict[seq_id]

        # add SeqRecord annotations
        current_SeqRecord_annotations = current_SeqRecord.annotations
        current_SeqRecord_annotations['date'] =                 (datetime.now().strftime('%d-%b-%Y')).upper()
        current_SeqRecord_annotations['accession'] =            ''
        current_SeqRecord_annotations['version'] =              ''
        current_SeqRecord_annotations['keywords'] =             ['.']
        current_SeqRecord_annotations['source'] =               prefix
        current_SeqRecord_annotations['organism'] =             prefix
        current_SeqRecord_annotations['taxonomy'] =             ['Unclassified']
        current_SeqRecord_annotations['comment'] =              '.'

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
                sequence_nc = str(Seq(id_to_sequence_dict[seq_id][cds_start-1:cds_end]).reverse_complement())

            # translate to aa sequence
            sequence_aa = str(SeqRecord(Seq(sequence_nc)).seq.translate(table=transl_table))

            # remove * at the end
            sequence_aa = sequence_aa[:-1]

            # export nc and aa sequences
            #export_dna_record(sequence_nc, locus_tag_id, '', bin_ffn_file_handle)
            bin_ffn_file_handle.write('>%s\n' % locus_tag_id)
            bin_ffn_file_handle.write('%s\n' % sequence_nc)
            #export_aa_record(sequence_aa, locus_tag_id, '', bin_faa_file_handle)
            bin_faa_file_handle.write('>%s\n' % locus_tag_id)
            bin_faa_file_handle.write('%s\n' % sequence_aa)

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

        # export to gbk file
        SeqIO.write(current_SeqRecord, bin_gbk_file_handle, 'genbank')

    bin_gbk_file_handle.close()
    bin_ffn_file_handle.close()
    bin_faa_file_handle.close()


def prodigal_worker(argument_list):

    input_genome                = argument_list[0]
    input_genome_folder         = argument_list[1]
    pwd_prodigal_exe            = argument_list[2]
    nonmeta_mode                = argument_list[3]
    pwd_prodigal_output_folder  = argument_list[4]

    # prepare command (according to Prokka)
    input_genome_basename, input_genome_ext = os.path.splitext(input_genome)
    pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    pwd_output_sco = '%s/%s.sco' % (pwd_prodigal_output_folder, input_genome_basename)

    prodigal_cmd_meta = '%s -f sco -q -c -m -g 11 -p meta -i %s -o %s' % (pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)
    prodigal_cmd_nonmeta = '%s -f sco -q -c -m -g 11 -i %s -o %s' % (pwd_prodigal_exe, pwd_input_genome, pwd_output_sco)

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


def remove_empty_element(list_in):

    list_out = []
    for each_element in list_in:
        if each_element != '':
            list_out.append(each_element)

    return list_out


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
            gap_percent = (dash_number / float(sequence_number)) * 100

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
                each_aa_percent = (current_column.count(each_aa) / float(sequence_number)) * 100
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


def parallel_blastn_worker(argument_list):
    query_file = argument_list[0]
    pwd_query_folder = argument_list[1]
    pwd_blast_db = argument_list[2]
    pwd_blast_result_folder = argument_list[3]
    blast_parameters = argument_list[4]
    pwd_blastn_exe = argument_list[5]

    pwd_blast_result_file = '%s/%s_blastn.tab' % (pwd_blast_result_folder, '.'.join(query_file.split('.')[:-1]))
    blastn_cmd = '%s -query %s/%s -db %s -out %s %s' % (pwd_blastn_exe,
                                                        pwd_query_folder,
                                                        query_file,
                                                        pwd_blast_db,
                                                        pwd_blast_result_file,
                                                        blast_parameters)
    os.system(blastn_cmd)


def PI(args, config_dict):

    # read in arguments
    input_genome_folder =   args['i']
    GTDB_output_file =      args['taxon']
    output_folder =         args['o']
    output_prefix =         args['p']
    grouping_levels =       args['r']
    grouping_file =         args['g']
    file_extension =        args['x']
    nonmeta_mode =          args['nonmeta']
    num_threads =           args['t']
    keep_quiet =            args['quiet']
    force_overwrite =       args['force']
    noblast =               args['noblast']

    # read in config file
    path_to_hmm =           config_dict['path_to_hmm']
    pwd_makeblastdb_exe =   config_dict['makeblastdb']
    pwd_blastn_exe =        config_dict['blastn']
    pwd_prodigal_exe =      config_dict['prodigal']
    pwd_hmmsearch_exe =     config_dict['hmmsearch']
    pwd_hmmfetch_exe =      config_dict['hmmfetch']
    pwd_hmmalign_exe =      config_dict['hmmalign']
    pwd_hmmstat_exe =       config_dict['hmmstat']
    pwd_fasttree_exe =      config_dict['fasttree']

    warnings.filterwarnings("ignore")

    minimal_cov_in_msa          = 50
    min_consensus_in_msa        = 25
    rank_abbre_dict             = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
    rank_abbre_dict_plural      = {'d': 'domains', 'p': 'phyla', 'c': 'classes', 'o': 'orders', 'f': 'families', 'g': 'genera', 's': 'species'}
    rank_to_position_dict       = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
    blast_parameters            = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % 1


    ######################################## check input file and dependencies #########################################

    if input_genome_folder[-1] == '/':
        input_genome_folder = input_genome_folder[:-1]

    # check whether executables exist
    program_list = [pwd_makeblastdb_exe, pwd_blastn_exe, pwd_prodigal_exe, pwd_hmmsearch_exe, pwd_hmmfetch_exe, pwd_hmmalign_exe, pwd_hmmstat_exe, pwd_fasttree_exe]
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)
    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


    # check parameters
    if (grouping_levels is not None) and (GTDB_output_file is None):
        print('Taxonomic classifications not detected, program exited')
        exit()

    # check whether input genome exist
    input_genome_file_re = '%s/*.%s' % (input_genome_folder, file_extension)
    input_genome_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_file_re)]
    input_genome_basename_list = ['.'.join(i.split('.')[:-1]) for i in input_genome_file_name_list]

    if input_genome_file_name_list == []:
        print('No input genome detected, program exited!')
        exit()


    # check the length of sequence id, exit if longer than 22bp
    # long_seq_id_genomes_list = []
    # for input_genome in input_genome_file_name_list:
    #     pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    #     for seq_record in SeqIO.parse(pwd_input_genome, 'fasta'):
    #         if len(seq_record.id) > 22:
    #             if input_genome not in long_seq_id_genomes_list:
    #                 long_seq_id_genomes_list.append(input_genome)
    #
    # if long_seq_id_genomes_list != []:
    #     print('Sequence id in the following genomes are longer than 22 letters, please shorten them!')
    #     for long_seq_id_genome in sorted(long_seq_id_genomes_list):
    #         print(long_seq_id_genome)
    #     print('Sequence id in the above genomes are longer than 22 letters, please shorten them!')
    #     exit()


    if grouping_levels is None:
        grouping_levels = 'x'

    # define folder name
    MetaCHIP_wd =    '%s_MetaCHIP_wd'        % (output_prefix)
    if output_folder is not None:
        MetaCHIP_wd = output_folder

    pwd_log_folder = '%s/%s_%s_log_files'    % (MetaCHIP_wd, output_prefix, grouping_levels)
    pwd_log_file =   '%s/%s_%s_PI_%s.log'    % (pwd_log_folder, output_prefix, grouping_levels, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))
    pwd_ignored_taxonomic_rank_file = '%s/ignored_taxonomic_rank.txt' % MetaCHIP_wd

    if (os.path.isdir(MetaCHIP_wd) is True) and (force_overwrite is False):
        print('MetaCHIP working directory detected, program exited!')
        exit()
    else:
        force_create_folder(MetaCHIP_wd)
        force_create_folder(pwd_log_folder)


    ############################################ read GTDB output into dict  ###########################################

    genomes_with_grouping = set()
    taxon_2_genome_dict_of_dict = {}
    ignored_genome_num = 0
    if grouping_levels == 'x':

        # read in grouping file
        group_id_2_genome_dict = {}
        for each_genome in open(grouping_file):
            each_genome_split = each_genome.strip().split(',')
            group_id = each_genome_split[0]
            genome_name = each_genome_split[1]

            genomes_with_grouping.add(genome_name)

            if group_id not in group_id_2_genome_dict:
                group_id_2_genome_dict[group_id] = [genome_name]
            else:
                group_id_2_genome_dict[group_id].append(genome_name)

        taxon_2_genome_dict_of_dict['x'] = group_id_2_genome_dict

    else:
        # read GTDB output into dict
        taxon_assignment_dict = {}
        for each_genome in open(GTDB_output_file):
            if not each_genome.startswith('user_genome'):
                each_split = each_genome.strip().split('\t')

                if len(each_split) == 1:
                    print('Unrecognisable %s, please make sure columns are tab separated, program exited!' % GTDB_output_file)
                    exit()

                bin_name = each_split[0]

                if bin_name in input_genome_basename_list:

                    assignment_full = []
                    if len(each_split) == 1:
                        assignment_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
                    elif (len(each_split) > 1) and (';' in each_split[1]):
                        assignment = each_split[1].split(';')
                        if len(assignment) == 7:
                            assignment_full = assignment
                        if len(assignment) == 6:
                            assignment_full = assignment + ['s__']
                        if len(assignment) == 5:
                            assignment_full = assignment + ['g__', 's__']
                        if len(assignment) == 4:
                            assignment_full = assignment + ['f__', 'g__', 's__']
                        if len(assignment) == 3:
                            assignment_full = assignment + ['o__', 'f__', 'g__', 's__']
                        if len(assignment) == 2:
                            assignment_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

                    elif (len(each_split) > 1) and (';' not in each_split[1]):
                        assignment_full = [each_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

                    # store in dict
                    taxon_assignment_dict[bin_name] = assignment_full

        # get all identified taxon at defined ranks
        for grouping_level in grouping_levels:

            taxon_2_genome_dict = {}

            specified_rank_pos = rank_to_position_dict[grouping_level]
            identified_taxon_list = []
            for each_TaxonAssign in taxon_assignment_dict:
                specified_rank_id = taxon_assignment_dict[each_TaxonAssign][specified_rank_pos]
                if specified_rank_id not in identified_taxon_list:
                    identified_taxon_list.append(specified_rank_id)

            # get the id of genomes assigned to each taxon at specified level
            for each_taxon in identified_taxon_list:
                genome_list = []
                for genome in taxon_assignment_dict:
                    if taxon_assignment_dict[genome][specified_rank_pos] == each_taxon:
                        genome_list.append(genome)
                taxon_2_genome_dict[each_taxon] = genome_list

            # get the number of ignored genome
            unclassified_symbol = '%s__' % grouping_level
            if unclassified_symbol in taxon_2_genome_dict:
                ignored_genome_num = len(taxon_2_genome_dict[unclassified_symbol])

            # report group number
            taxon_2_genome_dict.pop(unclassified_symbol, None)
            group_num = len(taxon_2_genome_dict)

            # for report and log
            sleep(0.5)
            report_and_log(('Input genomes grouped into %s %s.' % (group_num, rank_abbre_dict_plural[grouping_level])), pwd_log_file, keep_quiet)

            # report ignored genomes
            if ignored_genome_num > 0:
                sleep(0.5)
                report_and_log(('Ignored %s genome(s) for %s level HGT detection (unknown %s assignment).' % (ignored_genome_num, rank_abbre_dict[grouping_level], rank_abbre_dict[grouping_level])), pwd_log_file, keep_quiet)

            taxon_2_genome_dict_of_dict[grouping_level] = taxon_2_genome_dict

    taxon_2_genome_dict_of_dict_qualified = {}
    ignored_rank_list = []
    for each_rank in taxon_2_genome_dict_of_dict:
        each_rank_group_num = len(taxon_2_genome_dict_of_dict[each_rank])

        if each_rank_group_num > 1:
            taxon_2_genome_dict_of_dict_qualified[each_rank] = taxon_2_genome_dict_of_dict[each_rank]
        else:
            ignored_rank_list.append(each_rank)

    # write out ignored taxonomic rank
    if len(ignored_rank_list) > 0:
        pwd_ignored_taxonomic_rank_file_handle = open(pwd_ignored_taxonomic_rank_file, 'w')
        pwd_ignored_taxonomic_rank_file_handle.write('%s\n' % '\n'.join(ignored_rank_list))
        pwd_ignored_taxonomic_rank_file_handle.close()

    if len(ignored_rank_list) == len(taxon_2_genome_dict_of_dict):
        sleep(0.5)
        report_and_log(('Input genomes come from the same taxonomic group at all specified levels, program exited!'),pwd_log_file, keep_quiet)
        report_and_log(('Please note that file extension (e.g. fa, fasta) of the input genomes should NOT be included in the taxonomy or grouping file.'),pwd_log_file, keep_quiet)
        exit()
    else:
        for ignored_rank in ignored_rank_list:
            sleep(0.5)
            report_and_log(('Input genomes come from the same %s, ignored %s level HGT detection.' % (rank_abbre_dict[ignored_rank], rank_abbre_dict[ignored_rank])),pwd_log_file, keep_quiet)

    genome_for_HGT_detection_list = []
    for each_qualified_rank in taxon_2_genome_dict_of_dict_qualified:
        current_rank_group_to_genome_dict = taxon_2_genome_dict_of_dict_qualified[each_qualified_rank]
        for each_group in current_rank_group_to_genome_dict:
            genome_members = current_rank_group_to_genome_dict[each_group]
            for each_genome in genome_members:
                if each_genome not in genome_for_HGT_detection_list:
                    genome_for_HGT_detection_list.append('%s' % each_genome)

    sleep(0.5)
    report_and_log(('Total number of qualified genomes for HGT detection: %s.' % len(genome_for_HGT_detection_list)), pwd_log_file, keep_quiet)


    ############################################# define file/folder names #############################################

    combined_ffn_file =                  '%s_%s_combined_ffn.fasta'            % (output_prefix, grouping_levels)
    combined_faa_file =                  '%s_%s_combined_faa.fasta'            % (output_prefix, grouping_levels)
    prodigal_output_folder =             '%s_%s_prodigal_output'               % (output_prefix, grouping_levels)
    newick_tree_file =                   '%s_%s_SCG_tree.newick'               % (output_prefix, grouping_levels)
    SCG_tree_wd =                        '%s_%s_get_SCG_tree_wd'               % (output_prefix, grouping_levels)
    combined_alignment_file_tmp =        '%s_%s_SCG_tree_tmp.aln'              % (output_prefix, grouping_levels)
    combined_alignment_file =            '%s_%s_SCG_tree_cov%s_css%s.aln'      % (output_prefix, grouping_levels, minimal_cov_in_msa, min_consensus_in_msa)
    hmm_profile_sep_folder =             '%s_%s_hmm_profile_fetched'           % (output_prefix, grouping_levels)
    blast_db_folder =                    '%s_%s_blastdb'                       % (output_prefix, grouping_levels)
    blast_cmd_file =                     '%s_%s_blastn_commands.txt'           % (output_prefix, grouping_levels)
    blast_result_folder =                '%s_%s_blastn_results'                % (output_prefix, grouping_levels)

    pwd_combined_ffn_file =              '%s/%s/%s'                            % (MetaCHIP_wd, blast_db_folder, combined_ffn_file)
    pwd_combined_faa_file =              '%s/%s'                               % (MetaCHIP_wd, combined_faa_file)
    pwd_prodigal_output_folder =         '%s/%s'                               % (MetaCHIP_wd, prodigal_output_folder)
    pwd_combined_alignment_file_tmp =    '%s/%s/%s'                            % (MetaCHIP_wd, SCG_tree_wd, combined_alignment_file_tmp)
    pwd_hmm_profile_sep_folder =         '%s/%s/%s'                            % (MetaCHIP_wd, SCG_tree_wd, hmm_profile_sep_folder)
    pwd_combined_alignment_file =        '%s/%s'                               % (MetaCHIP_wd, combined_alignment_file)
    pwd_SCG_tree_wd =                    '%s/%s'                               % (MetaCHIP_wd, SCG_tree_wd)
    pwd_newick_tree_file =               '%s/%s'                               % (MetaCHIP_wd, newick_tree_file)
    pwd_blast_db_folder =                '%s/%s'                               % (MetaCHIP_wd, blast_db_folder)
    pwd_blast_result_folder =            '%s/%s'                               % (MetaCHIP_wd, blast_result_folder)
    pwd_blast_cmd_file =                 '%s/%s'                               % (MetaCHIP_wd, blast_cmd_file)


    ################################################### get grouping ###################################################

    if grouping_levels == 'x':
        pwd_grouping_file = grouping_file
        genomes_in_grouping_file = [i.strip().split(',')[1] for i in open(pwd_grouping_file)]
        if sorted(input_genome_basename_list) != sorted(genomes_in_grouping_file):
            report_and_log(('Genome ids provided in %s do not match genome files in %s, program exited!' % (pwd_grouping_file, input_genome_folder)), pwd_log_file, keep_quiet)
            report_and_log(('Please note that file extension (e.g. fa, fasta) of the input genomes should NOT be included in the grouping file.'), pwd_log_file, keep_quiet)
            exit()

    else:
        group_index_list = get_group_index_list()
        for grouping_level in taxon_2_genome_dict_of_dict_qualified:
            group_num = len(taxon_2_genome_dict_of_dict_qualified[grouping_level])
            grouping_file_name = '%s_grouping_%s%s.txt' % (output_prefix, grouping_level, group_num)
            pwd_grouping_file = '%s/%s' % (MetaCHIP_wd, grouping_file_name)

            grouping_file_handle = open(pwd_grouping_file, 'w')
            n = 0
            for each_taxon in taxon_2_genome_dict_of_dict_qualified[grouping_level]:
                group_id = group_index_list[n]
                for genome in taxon_2_genome_dict_of_dict_qualified[grouping_level][each_taxon]:
                    genomes_with_grouping.add(genome)
                    for_write = '%s,%s,%s\n'  % (group_id, genome, each_taxon)
                    grouping_file_handle.write(for_write)
                n += 1
            grouping_file_handle.close()

            # for report and log
            sleep(0.5)
            report_and_log(('Grouping file exported to: %s.' % grouping_file_name), pwd_log_file, keep_quiet)


    ######################################## run prodigal with multiprocessing #########################################

    # for report and log
    report_and_log(('Running Prodigal for %s qualified genomes with %s cores (1-3 minutes per genome per core).' % (len(genome_for_HGT_detection_list), num_threads)), pwd_log_file, keep_quiet)

    # create prodigal output folder
    os.mkdir(pwd_prodigal_output_folder)

    # prepare arguments for prodigal_worker
    list_for_multiple_arguments_Prodigal = []
    for input_genome in genome_for_HGT_detection_list:
        input_genome_with_extension = '%s.%s' % (input_genome, file_extension)
        list_for_multiple_arguments_Prodigal.append([input_genome_with_extension, input_genome_folder, pwd_prodigal_exe, nonmeta_mode, pwd_prodigal_output_folder])

    # run prodigal with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(prodigal_worker, list_for_multiple_arguments_Prodigal)
    pool.close()
    pool.join()


    ########################################### get species tree (hmmsearch) ###########################################

    # create wd
    os.mkdir(pwd_SCG_tree_wd)

    # for report and log
    report_and_log(('Get SCG tree: running hmmsearch with %s cores.' % num_threads), pwd_log_file, keep_quiet)

    # prepare arguments for hmmsearch_worker
    list_for_multiple_arguments_hmmsearch = []
    for faa_file_basename in genome_for_HGT_detection_list:
        list_for_multiple_arguments_hmmsearch.append([faa_file_basename, pwd_SCG_tree_wd, pwd_hmmsearch_exe, path_to_hmm, pwd_prodigal_output_folder])

    # run hmmsearch with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(hmmsearch_worker, list_for_multiple_arguments_hmmsearch)
    pool.close()
    pool.join()


    ############################################# get species tree (hmmalign) #############################################

    # for report and log
    report_and_log(('Get SCG tree: running hmmalign with %s cores.' % num_threads), pwd_log_file, keep_quiet)

    # fetch combined hmm profiles
    os.mkdir(pwd_hmm_profile_sep_folder)
    sep_combined_hmm(path_to_hmm, pwd_hmm_profile_sep_folder, pwd_hmmfetch_exe, pwd_hmmstat_exe)

    # Call hmmalign to align all single fasta files with hmms
    files = os.listdir(pwd_SCG_tree_wd)
    fastaFiles = [i for i in files if i.endswith('.fasta')]

    # prepare arguments for hmmalign_worker
    list_for_multiple_arguments_hmmalign = []
    for fastaFile in fastaFiles:

        fastaFiles_basename = '.'.join(fastaFile.split('.')[:-1])
        list_for_multiple_arguments_hmmalign.append([fastaFiles_basename, pwd_SCG_tree_wd, pwd_hmm_profile_sep_folder, pwd_hmmalign_exe])

    # run hmmalign with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(hmmalign_worker, list_for_multiple_arguments_hmmalign)
    pool.close()
    pool.join()


    ################################### get species tree (Concatenating alignments) ####################################

    # for report and log
    report_and_log('Get SCG tree: concatenating alignments.', pwd_log_file, keep_quiet)

    # concatenating the single alignments
    concatAlignment = {}
    for element in genome_for_HGT_detection_list:
        concatAlignment[element] = ''

    # Reading all single alignment files and append them to the concatenated alignment
    files = os.listdir(pwd_SCG_tree_wd)
    fastaFiles = [i for i in files if i.endswith('.fasta')]
    for faa_file_basename in fastaFiles:
        fastaFile = pwd_SCG_tree_wd + '/' + faa_file_basename
        proteinSequence = {}
        alignmentLength = 0
        for seq_record_2 in SeqIO.parse(fastaFile, 'fasta'):
            proteinName = seq_record_2.id
            proteinSequence[proteinName] = str(seq_record_2.seq)
            alignmentLength = len(proteinSequence[proteinName])

        for element in genome_for_HGT_detection_list:
            if element in proteinSequence.keys():
                concatAlignment[element] += proteinSequence[element]
            else:
                concatAlignment[element] += '-' * alignmentLength

    # writing alignment to file
    file_out = open(pwd_combined_alignment_file_tmp, 'w')
    for element in genome_for_HGT_detection_list:
        file_out.write('>' + element + '\n' + concatAlignment[element] + '\n')
    file_out.close()

    # remove columns with low coverage and low consensus
    report_and_log(('Get SCG tree: removing columns from concatenated alignment represented by <%s%s of genomes and with amino acid consensus <%s%s.' % (minimal_cov_in_msa, '%', min_consensus_in_msa, '%')), pwd_log_file, keep_quiet)
    remove_low_cov_and_consensus_columns(pwd_combined_alignment_file_tmp, minimal_cov_in_msa, min_consensus_in_msa, pwd_combined_alignment_file)


    ########################################### get species tree (fasttree) ############################################

    # for report and log
    report_and_log('Get SCG tree: running FastTree.', pwd_log_file, keep_quiet)

    # calling fasttree for tree calculation
    fasttree_cmd = '%s -quiet %s > %s 2>/dev/null' % (pwd_fasttree_exe, pwd_combined_alignment_file, pwd_newick_tree_file)
    os.system(fasttree_cmd)

    # for report and log
    report_and_log(('SCG tree exported to: %s.' % newick_tree_file), pwd_log_file, keep_quiet)


    ############################################### run all vs all blastn ##############################################

    # get combined faa file
    os.system('cat %s/*.faa > %s' % (pwd_prodigal_output_folder, pwd_combined_faa_file))

    # create folder
    os.mkdir(pwd_blast_db_folder)
    os.mkdir(pwd_blast_result_folder)

    # run makeblastdb
    report_and_log(('Making blast database.'), pwd_log_file, keep_quiet)
    os.system('cat %s/*.ffn > %s/%s' % (pwd_prodigal_output_folder, pwd_blast_db_folder, combined_ffn_file))

    makeblastdb_cmd = '%s -in %s/%s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_blast_db_folder, combined_ffn_file)
    os.system(makeblastdb_cmd)

    # check db files

    # prepare arguments list for parallel_blastn_worker
    ffn_file_list = ['%s.ffn' % i for i in genome_for_HGT_detection_list]

    pwd_blast_cmd_file_handle = open(pwd_blast_cmd_file, 'w')
    list_for_multiple_arguments_blastn = []
    for ffn_file in ffn_file_list:
        list_for_multiple_arguments_blastn.append([ffn_file, pwd_prodigal_output_folder, pwd_combined_ffn_file, pwd_blast_result_folder, blast_parameters, pwd_blastn_exe])
        blastn_cmd = '%s -query %s/%s -db %s -out %s/%s %s' % (pwd_blastn_exe, pwd_prodigal_output_folder, ffn_file, pwd_combined_ffn_file, pwd_blast_result_folder, '%s_blastn.tab' % '.'.join(ffn_file.split('.')[:-1]), blast_parameters)
        pwd_blast_cmd_file_handle.write('%s\n' % blastn_cmd)
    pwd_blast_cmd_file_handle.close()

    report_and_log(('Blastn commands exported to: %s.' % blast_cmd_file), pwd_log_file, keep_quiet)

    if noblast is False:

        report_and_log(('Running blastn for %s qualified genomes with %s cores.' % (len(genome_for_HGT_detection_list), num_threads)), pwd_log_file, keep_quiet)

        # run blastn with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(parallel_blastn_worker, list_for_multiple_arguments_blastn)
        pool.close()
        pool.join()

        report_and_log(('Blast results exported to: %s.' % pwd_blast_result_folder), pwd_log_file, keep_quiet)

    if noblast is False:
        report_and_log('PI step done!', pwd_log_file, keep_quiet)
    else:
        report_and_log('PI step done!', pwd_log_file, keep_quiet)
        report_and_log(('All-vs-all blastn disabled, please run blastn with commands in: %s before the BP module.' % blast_cmd_file), pwd_log_file, keep_quiet)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # arguments for PI
    parser.add_argument('-i',       required=True,                       help='input genome folder')
    parser.add_argument('-taxon',   required=False,                      help='taxonomic classification of input genomes')
    parser.add_argument('-o',       required=False, default=None,        help='output folder (default: current working directory)')
    parser.add_argument('-p',       required=True,                       help='output prefix')
    parser.add_argument('-r',       required=False, default=None,        help='grouping rank, choose from p (phylum), c (class), o (order), f (family), g (genus) or any combination of them')
    parser.add_argument('-g',       required=False, default=None,        help='grouping file')
    parser.add_argument('-x',       required=False, default='fasta',     help='file extension')
    parser.add_argument('-nonmeta', required=False, action="store_true", help='provide if input genomes are NOT metagenome-assembled genomes')
    parser.add_argument('-t',       required=False, type=int, default=1, help='number of threads, default: 1')
    parser.add_argument('-quiet',   required=False, action="store_true", help='not report progress')
    parser.add_argument('-force',   required=False, action="store_true", help='force overwrite existing results')
    parser.add_argument('-noblast', required=False, action="store_true", help='skip running all-vs-all blastn, provide if you have other ways (e.g. with job scripts) to speed up the blastn step')

    args = vars(parser.parse_args())

    PI(args, config_dict)

