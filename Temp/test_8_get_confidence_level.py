#!/usr/bin/env python
import os
import math
import matplotlib as mpl
mpl.use('Agg')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord



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


def get_confidence_level(folder_name, flanking_length, calculation_step, pwd_blastn_exe):

    # get the number of steps
    step_number = int(math.ceil(flanking_length / float(calculation_step)))
    print("step_number:")
    print(step_number)
    # define file name
    recipient_gene = folder_name.split('___')[0]
    donor_gene = folder_name.split('___')[1]
    file_recipient_gene_3000_gbk = '%s_%sbp.gbk' % (recipient_gene, flanking_length)
    file_donor_gene_3000_gbk = '%s_%sbp.gbk' % (donor_gene, flanking_length)

    # read in recipient contig
    recipient_contig_record = SeqIO.read(file_recipient_gene_3000_gbk, 'genbank')
    recipient_contig_seq = recipient_contig_record.seq
    print('\nrecipient_gene_location:')
    #print(recipient_contig_seq)
    gene_r_start = 0
    gene_r_end = 0
    gene_r_strand = 0
    for gene_r in recipient_contig_record.features:
        if 'locus_tag' in gene_r.qualifiers:
            if gene_r.qualifiers['locus_tag'][0] == recipient_gene:
                gene_r_start = int(gene_r.location.start)
                gene_r_end = int(gene_r.location.end)
                gene_r_strand = int(gene_r.location.strand)
                print(gene_r_start)
                print(gene_r_end)
                print(gene_r_strand)
    print('recipient_contig_length: %s' % len(recipient_contig_seq))

    # read in donor contig
    donor_contig_record = SeqIO.read(file_donor_gene_3000_gbk, 'genbank')
    donor_contig_seq = donor_contig_record.seq
    print('\ndonor_gene_location:')
    gene_d_start = 0
    gene_d_end = 0
    gene_d_strand = 0
    for gene_d in donor_contig_record.features:
        if 'locus_tag' in gene_d.qualifiers:
            if gene_d.qualifiers['locus_tag'][0] == donor_gene:
                gene_d_start = int(gene_d.location.start)
                gene_d_end = int(gene_d.location.end)
                gene_d_strand = int(gene_d.location.strand)
                print(gene_d_start)
                print(gene_d_end)
                print(gene_d_strand)
    print('donor_contig_length: %s' % len(donor_contig_seq))
    # get the length of subsequences
    subseq_length_dict = {}

    # get subsequences of recipient contig left flanking region
    rlf_subs_list = []
    n_rlf = 1
    while n_rlf <= step_number:
        current_bp = gene_r_start - n_rlf * calculation_step
        if current_bp < 0:
            current_bp = 0
        rlf_subs_list.append(current_bp)
        n_rlf += 1
    print('\nrlf_subs_list:')
    print(rlf_subs_list)
    n_rlf = 0
    while n_rlf <= step_number-1:
        if n_rlf == 0:
            current_rlf_seq = recipient_contig_seq[rlf_subs_list[n_rlf]:gene_r_start]
        else:
            current_rlf_seq = recipient_contig_seq[rlf_subs_list[n_rlf]:rlf_subs_list[n_rlf-1]]

        current_rlf_id = '%s_rlf%s' % (recipient_gene, n_rlf + 1)
        if len(current_rlf_seq) > 0:
            current_rlf_description = ''
            current_rlf_handle = '%s/%s.fasta' % (os.getcwd(), current_rlf_id)
            export_dna_record(current_rlf_seq, current_rlf_id, current_rlf_description, current_rlf_handle)
        subseq_length_dict[current_rlf_id] = len(current_rlf_seq)
        print(current_rlf_id)
        #print(current_rlf_seq)
        print(len(current_rlf_seq))
        n_rlf += 1

    # get subsequences of recipient contig right flanking region
    rrf_subs_list = []
    n_rrf = 1
    while n_rrf <= step_number:
        current_bp = gene_r_end + n_rrf * calculation_step
        if current_bp > len(recipient_contig_seq):
            current_bp = len(recipient_contig_seq)
        rrf_subs_list.append(current_bp)
        n_rrf += 1
    print('\nrrf_subs_list')
    print(rrf_subs_list)
    n_rrf = 0
    while n_rrf <= step_number-1:
        if n_rrf == 0:
            current_rrf_seq = recipient_contig_seq[gene_r_end:rrf_subs_list[n_rrf]]
        else:
            current_rrf_seq = recipient_contig_seq[rrf_subs_list[n_rrf-1]:rrf_subs_list[n_rrf]]
        current_rrf_id = '%s_rrf%s' % (recipient_gene, n_rrf + 1)
        if len(current_rrf_seq) > 0:
            current_rrf_description = ''
            current_rrf_handle = '%s/%s.fasta' % (os.getcwd(), current_rrf_id)
            export_dna_record(current_rrf_seq, current_rrf_id, current_rrf_description, current_rrf_handle)
        subseq_length_dict[current_rrf_id] = len(current_rrf_seq)
        print(current_rrf_id)
        print(len(current_rrf_seq))
        #print(current_rrf_seq)
        n_rrf += 1

    # get subsequences of donor contig left flanking region
    dlf_subs_list = []
    n_dlf = 1
    while n_dlf <= step_number:
        current_bp = gene_d_start - n_dlf * calculation_step
        if current_bp < 0:
            current_bp = 0
        dlf_subs_list.append(current_bp)
        n_dlf += 1
    print('\ndlf_subs_list:')
    print(dlf_subs_list)
    n_dlf = 0
    while n_dlf <= step_number-1:
        if n_dlf == 0:
            current_dlf_seq = donor_contig_seq[dlf_subs_list[n_dlf]:gene_d_start]
        else:
            current_dlf_seq = donor_contig_seq[dlf_subs_list[n_dlf]:dlf_subs_list[n_dlf-1]]
        current_dlf_id = '%s_dlf%s' % (donor_gene, n_dlf + 1)
        if len(current_dlf_seq) >0:
            current_dlf_description = ''
            current_dlf_handle = '%s/%s.fasta' % (os.getcwd(), current_dlf_id)
            export_dna_record(current_dlf_seq, current_dlf_id, current_dlf_description, current_dlf_handle)
        subseq_length_dict[current_dlf_id] = len(current_dlf_seq)
        print(current_dlf_id)
        print(len(current_dlf_seq))
        n_dlf += 1

    # get subsequences of donor contig right flanking region
    drf_subs_list = []
    n_drf = 1
    while n_drf <= step_number:
        current_bp = gene_d_end + n_drf * calculation_step
        if current_bp > len(donor_contig_seq):
            current_bp = len(donor_contig_seq)
        drf_subs_list.append(current_bp)
        n_drf += 1
    print('\ndrf_subs_list:')
    print(drf_subs_list)
    n_drf = 0
    while n_drf <= step_number-1:
        if n_drf == 0:
            current_drf_seq = donor_contig_seq[gene_d_end:drf_subs_list[n_drf]]
        else:
            current_drf_seq = donor_contig_seq[drf_subs_list[n_drf-1]:drf_subs_list[n_drf]]
        current_drf_id = '%s_drf%s' % (donor_gene, n_drf + 1)
        if len(current_drf_seq) > 0:
            current_drf_description = ''
            current_drf_handle = '%s/%s.fasta' % (os.getcwd(), current_drf_id)
            export_dna_record(current_drf_seq, current_drf_id, current_drf_description, current_drf_handle)
        subseq_length_dict[current_drf_id] = len(current_drf_seq)
        print(current_drf_id)
        print(len(current_drf_seq))
        n_drf += 1
    print('')
    print(subseq_length_dict)

    # run pair-wise blast
    blast_parameters = '-evalue 1e-5 -outfmt 6 -task blastn'
    n = 1
    alignment_length_cutoff = round(float(calculation_step) / 10)
    recipient_left_iden_profile = []
    recipient_right_iden_profile = []
    while n <= step_number:
        query_l = '%s/%s_rlf%s.fasta' % (os.getcwd(), recipient_gene, n)
        query_l_length = subseq_length_dict['%s_rlf%s' % (recipient_gene, n)]
        query_r = '%s/%s_rrf%s.fasta' % (os.getcwd(), recipient_gene, n)
        query_r_length = subseq_length_dict['%s_rrf%s' % (recipient_gene, n)]
        donor_l = '%s/%s_dlf%s.fasta' % (os.getcwd(), donor_gene, n)
        donor_l_length = subseq_length_dict['%s_dlf%s' % (donor_gene, n)]
        donor_r = '%s/%s_drf%s.fasta' % (os.getcwd(), donor_gene, n)
        donor_r_length = subseq_length_dict['%s_drf%s' % (donor_gene, n)]
        print()
        print('%s\t%s' % (query_l, query_l_length))
        print('%s\t%s' % (query_r, query_r_length))
        print('%s\t%s' % (donor_l, donor_l_length))
        print('%s\t%s' % (donor_r, donor_r_length))

        # run pair-wise blast
        if gene_r_strand == gene_d_strand:
            output_1 = '%s/%s_rlf%s___%s_dlf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            output_2 = '%s/%s_rrf%s___%s_drf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            command_blast_1 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_l, donor_l, output_1, blast_parameters)
            if (query_l_length > 0) and (donor_l_length > 0):
                os.system(command_blast_1)
            command_blast_2 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_r, donor_r, output_2, blast_parameters)
            if (query_r_length > 0) and (donor_r_length > 0):
                os.system(command_blast_2)

        if gene_r_strand != gene_d_strand:
            output_1 = '%s/%s_rlf%s___%s_drf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            output_2 = '%s/%s_rrf%s___%s_dlf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            command_blast_1 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_l, donor_r, output_1, blast_parameters)
            if (query_l_length > 0) and (donor_r_length > 0):
                os.system(command_blast_1)
            command_blast_2 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_r, donor_l, output_2, blast_parameters)
            if (query_r_length > 0) and (donor_l_length > 0):
                os.system(command_blast_2)

        # get recipient_left_iden_profile
        if not os.path.isfile(output_1):
            recipient_left_iden_profile.append(['No_sequence'])
        elif os.path.isfile(output_1) and (os.path.getsize(output_1) == 0):
            recipient_left_iden_profile.append(['No_similarity'])
        elif os.path.getsize(output_1) > 0:
            identity_l_list = []
            for each_line_1 in open(output_1):
                each_line_1_split = each_line_1.strip().split('\t')
                identity_1 = round(float(each_line_1_split[2]))
                alignment_length_1 = int(each_line_1_split[3])
                if alignment_length_1 >= alignment_length_cutoff:
                    identity_l_list.append(identity_1)
            if identity_l_list == []:
                recipient_left_iden_profile.append(['No_similarity'])
            else:
                recipient_left_iden_profile.append(identity_l_list)

        # get recipient_right_iden_profile
        if not os.path.isfile(output_2):
            recipient_right_iden_profile.append(['No_sequence'])
        elif os.path.isfile(output_2) and (os.path.getsize(output_2) == 0):
            recipient_right_iden_profile.append(['No_similarity'])
        elif os.path.getsize(output_2) > 0:
            identity_r_list = []
            for each_line_2 in open(output_2):
                each_line_2_split = each_line_2.strip().split('\t')
                identity_2 = round(float(each_line_2_split[2]))
                alignment_length_2 = int(each_line_2_split[3])
                if alignment_length_2 >= alignment_length_cutoff:
                    identity_r_list.append(identity_2)
            if identity_r_list == []:
                recipient_right_iden_profile.append(['No_similarity'])
            else:
                recipient_right_iden_profile.append(identity_r_list)
        n += 1

    print('recipient_left_iden_profile: %s' % recipient_left_iden_profile)
    print('recipient_right_iden_profile: %s' % recipient_right_iden_profile)

    # analyze flanking profile
    recipient_left_iden_profile_uniq = uniq_list(recipient_left_iden_profile)
    recipient_right_iden_profile_uniq = uniq_list(recipient_right_iden_profile)
    print(recipient_left_iden_profile_uniq)
    print(recipient_right_iden_profile_uniq)
    confidence_level = ''

    if (['No_sequence'] not in recipient_left_iden_profile) and (['No_sequence'] not in recipient_right_iden_profile):
        if (recipient_left_iden_profile[-1] == ['No_similarity']) and (recipient_right_iden_profile[-1] == ['No_similarity']): # non-end match
            confidence_level = 'high'
        elif (recipient_left_iden_profile[-1] != ['No_similarity']) or (recipient_right_iden_profile[-1] != ['No_similarity']): # matched to the end
            confidence_level = 'low'
        else:
            confidence_level = 'unknown'

    elif (['No_sequence'] in recipient_left_iden_profile) and (['No_sequence'] not in recipient_right_iden_profile):
        if (recipient_left_iden_profile_uniq[-2] == ['No_similarity']) and (recipient_right_iden_profile_uniq[-1] == ['No_similarity']):
            confidence_level = 'high'
        else:
            confidence_level = 'low'

    elif (['No_sequence'] not in recipient_left_iden_profile) and (['No_sequence'] in recipient_right_iden_profile):
        if (recipient_left_iden_profile_uniq[-1] == ['No_similarity']) and (recipient_right_iden_profile_uniq[-2] == ['No_similarity']):
            confidence_level = 'high'
        else:
            confidence_level = 'low'

    elif (['No_sequence'] in recipient_left_iden_profile) and (['No_sequence'] in recipient_right_iden_profile):
        if (recipient_left_iden_profile_uniq[-2] == ['No_similarity']) and (recipient_right_iden_profile_uniq[-2] == ['No_similarity']):
            confidence_level = 'high'
        else:
            confidence_level = 'low'

    return confidence_level


wd = '/Users/songweizhi/Desktop/MetaCHIP_wd/optimize_confidence_level'
candidate = 'AA_Refined_13_00139___PA_Refined_11_01369'
# candidate = 'AA_Refined_13_00597___DP_Refined_163_00555'
# candidate = 'AA_Refined_109_01081___AA_Refined_96_01072'
# candidate = 'AA_Refined_116_00628___PC_Refined_57_00060'
# candidate = 'AA_Refined_117_00923___HO_Refined_109_00124'

os.chdir('%s/%s' % (wd, candidate))
confidence_level = get_confidence_level(candidate, 3000, 300, 'blastn')
print(confidence_level)
