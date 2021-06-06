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
import warnings
warnings.filterwarnings("ignore")

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


seq_file      = '/Users/songweizhi/Desktop/gbk_test/SH_ER_070217_Refined_44.fasta'
sco_file      = '/Users/songweizhi/Desktop/gbk_test/SH_ER_070217_Refined_44.sco'
prefix        = 'SH_ER_070217_Refined_44'
output_folder = '/Users/songweizhi/Desktop/gbk_test'
prodigal_parser(seq_file, sco_file, prefix, output_folder)


'''
prodigal -f sco -q -c -m -g 11 -p meta -i SH_ER_070217_Refined_44.fasta -o SH_ER_070217_Refined_44.sco
'''
