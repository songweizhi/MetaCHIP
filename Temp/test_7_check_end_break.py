#!/usr/bin/env python

import os
import math
import shutil
import argparse
import configparser
from sys import stdout
from time import sleep
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm


def uniq_list(input_list):
    output_list = []
    for each_element in input_list:
        if each_element not in output_list:
            output_list.append(each_element)
    return output_list


def export_dna_record(gene_seq, gene_id, gene_description, pwd_output_file):
    output_handle = open(pwd_output_file, 'w')
    seq_object = Seq(str(gene_seq), IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')
    output_handle.close()


def check_end_break(folder_name, flanking_length, end_seq_length, pwd_blastn_exe):
    # define file name
    recipient_gene = folder_name.split('___')[0]
    donor_gene = folder_name.split('___')[1]
    file_recipient_gene_3000_gbk = '%s_%sbp.gbk' % (recipient_gene, flanking_length)
    file_donor_gene_3000_gbk = '%s_%sbp.gbk' % (donor_gene, flanking_length)

    # read in recipient/donor contig
    recipient_contig_record = SeqIO.read(file_recipient_gene_3000_gbk, 'genbank')
    recipient_contig_seq = recipient_contig_record.seq
    donor_contig_record = SeqIO.read(file_donor_gene_3000_gbk, 'genbank')
    donor_contig_seq = donor_contig_record.seq

    # get ending sequence of the recipient and donor genes
    ending_seq_description = ''

    # export recipient_left_end_seq
    recipient_left_end_seq = recipient_contig_seq[0:end_seq_length]
    recipient_left_end_id = '%s_le%s' % (recipient_gene, end_seq_length)
    recipient_left_end_handle = '%s/%s.fasta' % (os.getcwd(), recipient_left_end_id)
    export_dna_record(recipient_left_end_seq, recipient_left_end_id, ending_seq_description, recipient_left_end_handle)

    # export recipient_right_end_seq
    recipient_right_end_seq = recipient_contig_seq[len(recipient_contig_seq) - end_seq_length:]
    recipient_right_end_id = '%s_re%s' % (recipient_gene, end_seq_length)
    recipient_right_end_handle = '%s/%s.fasta' % (os.getcwd(), recipient_right_end_id)
    export_dna_record(recipient_right_end_seq, recipient_right_end_id, ending_seq_description,
                      recipient_right_end_handle)

    # export donor_left_end_seq
    donor_left_end_seq = donor_contig_seq[0:end_seq_length]
    donor_left_end_id = '%s_le%s' % (donor_gene, end_seq_length)
    donor_left_end_handle = '%s/%s.fasta' % (os.getcwd(), donor_left_end_id)
    export_dna_record(donor_left_end_seq, donor_left_end_id, ending_seq_description, donor_left_end_handle)

    # export donor_right_end_seq
    donor_right_end_seq = donor_contig_seq[len(donor_contig_seq) - end_seq_length:]
    donor_right_end_id = '%s_re%s' % (donor_gene, end_seq_length)
    donor_right_end_handle = '%s/%s.fasta' % (os.getcwd(), donor_right_end_id)
    export_dna_record(donor_right_end_seq, donor_right_end_id, ending_seq_description, donor_right_end_handle)

    # run blastn between ending sequences:
    blast_parameters = '-evalue 1e-5 -outfmt 6 -task blastn'
    output_rle_dle = '%s/%s_rle___%s_dle.tab' % (os.getcwd(), recipient_gene, donor_gene)
    output_rle_dre = '%s/%s_rle___%s_dre.tab' % (os.getcwd(), recipient_gene, donor_gene)
    output_rre_dle = '%s/%s_rre___%s_dle.tab' % (os.getcwd(), recipient_gene, donor_gene)
    output_rre_dre = '%s/%s_rre___%s_dre.tab' % (os.getcwd(), recipient_gene, donor_gene)
    compare_end_blast_rle_dle = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_left_end_handle, donor_left_end_handle, output_rle_dle, blast_parameters)
    compare_end_blast_rle_dre = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_left_end_handle, donor_right_end_handle, output_rle_dre, blast_parameters)
    compare_end_blast_rre_dle = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_right_end_handle, donor_left_end_handle, output_rre_dle, blast_parameters)
    compare_end_blast_rre_dre = '%s -query %s -subject %s -out %s %s' % (
    pwd_blastn_exe, recipient_right_end_handle, donor_right_end_handle, output_rre_dre, blast_parameters)
    os.system(compare_end_blast_rle_dle)
    os.system(compare_end_blast_rle_dre)
    os.system(compare_end_blast_rre_dle)
    os.system(compare_end_blast_rre_dre)

    match_files = [output_rle_dle, output_rle_dre, output_rre_dle, output_rre_dre]
    match_profile = []
    for each_match in match_files:
        match_dir_list = []
        for each_hit in open(each_match):
            each_hit_split = each_hit.strip().split('\t')
            qstart = int(each_hit_split[6])
            qend = int(each_hit_split[7])
            sstart = int(each_hit_split[8])
            send = int(each_hit_split[9])
            q_direction = qend - qstart
            s_direction = send - sstart
            q_dir = ''
            if q_direction > 0:
                q_dir = 'forward'
            if q_direction < 0:
                q_dir = 'backward'
            s_dir = ''
            if s_direction > 0:
                s_dir = 'forward'
            if s_direction < 0:
                s_dir = 'backward'
            match_dir = ''
            if q_dir == s_dir:
                match_dir = 'same direction'
            if q_dir != s_dir:
                match_dir = 'opposite direction'
            match_dir_list.append(match_dir)
        match_dir_list_uniq = []
        for each in match_dir_list:
            if each not in match_dir_list_uniq:
                match_dir_list_uniq.append(each)
        match_profile.append(match_dir_list_uniq)

    break_end = ''
    if (match_profile[0] == ['opposite direction']) or (match_profile[1] == ['same direction']) or (match_profile[2] == ['same direction']) or (match_profile[3] == ['opposite direction']):
        break_end = True
    else:
        break_end = False

    return break_end


folder_name = 'CF_Refined_15_00035___HO_Refined_79_00184'
end_break = check_end_break(folder_name, 3000, 1000, 'blastn')
print(end_break)











