import os
import math
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


def get_match_category(folder_name, flanking_length, calculation_step, pwd_blastn_exe):

    # get the number of steps
    step_number = int(math.ceil(flanking_length / float(calculation_step)))

    # define file name
    recipient_gene = folder_name.split('___')[0]
    donor_gene = folder_name.split('___')[1]
    file_recipient_gene_3000_gbk = '%s_%sbp.gbk' % (recipient_gene, flanking_length)
    file_donor_gene_3000_gbk = '%s_%sbp.gbk' % (donor_gene, flanking_length)

    # read in recipient contig
    recipient_contig_record = SeqIO.read(file_recipient_gene_3000_gbk, 'genbank')
    recipient_contig_seq = recipient_contig_record.seq
    gene_r_start = 0
    gene_r_end = 0
    gene_r_strand = 0
    for gene_r in recipient_contig_record.features:
        if 'locus_tag' in gene_r.qualifiers:
            if gene_r.qualifiers['locus_tag'][0] == recipient_gene:
                gene_r_start = int(gene_r.location.start)
                gene_r_end = int(gene_r.location.end)
                gene_r_strand = int(gene_r.location.strand)

    # read in donor contig
    donor_contig_record = SeqIO.read(file_donor_gene_3000_gbk, 'genbank')
    donor_contig_seq = donor_contig_record.seq
    gene_d_start = 0
    gene_d_end = 0
    gene_d_strand = 0
    for gene_d in donor_contig_record.features:
        if 'locus_tag' in gene_d.qualifiers:
            if gene_d.qualifiers['locus_tag'][0] == donor_gene:
                gene_d_start = int(gene_d.location.start)
                gene_d_end = int(gene_d.location.end)
                gene_d_strand = int(gene_d.location.strand)

    # get subsequences of recipient contig left flanking region
    rlf_subs_list = []
    n_rlf = 1
    while n_rlf <= step_number:
        current_bp = gene_r_start - n_rlf * calculation_step
        rlf_subs_list.append(current_bp)
        n_rlf += 1

    n_rlf = 0
    while n_rlf <= step_number-1:
        if n_rlf == 0:
            current_rlf_seq = recipient_contig_seq[rlf_subs_list[n_rlf]:gene_r_start]
        else:
            current_rlf_seq = recipient_contig_seq[rlf_subs_list[n_rlf]:rlf_subs_list[n_rlf-1]]
        current_rlf_id = '%s_rlf%s' % (recipient_gene, n_rlf + 1)
        current_rlf_description = ''
        current_rlf_handle = '%s/%s.fasta' % (os.getcwd(), current_rlf_id)
        export_dna_record(current_rlf_seq, current_rlf_id, current_rlf_description, current_rlf_handle)
        n_rlf += 1

    # get subsequences of recipient contig right flanking region
    rrf_subs_list = []
    n_rrf = 1
    while n_rrf <= step_number:
        current_bp = gene_r_end + n_rrf * calculation_step
        rrf_subs_list.append(current_bp)
        n_rrf += 1

    n_rrf = 0
    while n_rrf <= step_number-1:
        if n_rrf == 0:
            current_rrf_seq = recipient_contig_seq[gene_r_end:rrf_subs_list[n_rrf]]
        else:
            current_rrf_seq = recipient_contig_seq[rrf_subs_list[n_rrf-1]:rrf_subs_list[n_rrf]]
        current_rrf_id = '%s_rrf%s' % (recipient_gene, n_rrf + 1)
        current_rrf_description = ''
        current_rrf_handle = '%s/%s.fasta' % (os.getcwd(), current_rrf_id)
        export_dna_record(current_rrf_seq, current_rrf_id, current_rrf_description, current_rrf_handle)
        n_rrf += 1

    # get subsequences of donor contig left flanking region
    dlf_subs_list = []
    n_dlf = 1
    while n_dlf <= step_number:
        current_bp = gene_d_start - n_dlf * calculation_step
        dlf_subs_list.append(current_bp)
        n_dlf += 1

    n_dlf = 0
    while n_dlf <= step_number-1:
        if n_dlf == 0:
            current_dlf_seq = donor_contig_seq[dlf_subs_list[n_dlf]:gene_d_start]
        else:
            current_dlf_seq = donor_contig_seq[dlf_subs_list[n_dlf]:dlf_subs_list[n_dlf-1]]
        current_dlf_id = '%s_dlf%s' % (donor_gene, n_dlf + 1)
        current_dlf_description = ''
        current_dlf_handle = '%s/%s.fasta' % (os.getcwd(), current_dlf_id)
        export_dna_record(current_dlf_seq, current_dlf_id, current_dlf_description, current_dlf_handle)
        n_dlf += 1

    # get subsequences of donor contig right flanking region
    drf_subs_list = []
    n_drf = 1
    while n_drf <= step_number:
        current_bp = gene_d_end + n_drf * calculation_step
        drf_subs_list.append(current_bp)
        n_drf += 1

    n_drf = 0
    while n_drf <= step_number-1:
        if n_drf == 0:
            current_drf_seq = donor_contig_seq[gene_d_end:drf_subs_list[n_drf]]
        else:
            current_drf_seq = donor_contig_seq[drf_subs_list[n_drf-1]:drf_subs_list[n_drf]]
        current_drf_id = '%s_drf%s' % (donor_gene, n_drf + 1)
        current_drf_description = ''
        current_drf_handle = '%s/%s.fasta' % (os.getcwd(), current_drf_id)
        export_dna_record(current_drf_seq, current_drf_id, current_drf_description, current_drf_handle)
        n_drf += 1

    # run pair-wise blast
    blast_parameters = '-evalue 1e-5 -outfmt 6 -task blastn'
    n = 1
    alignment_length_cutoff = round(float(calculation_step) / 5)
    recipient_left_iden_profile = []
    recipient_right_iden_profile = []
    while n <= step_number:
        query_l = '%s/%s_rlf%s.fasta' % (os.getcwd(), recipient_gene, n)
        query_r = '%s/%s_rrf%s.fasta' % (os.getcwd(), recipient_gene, n)
        donor_l = '%s/%s_dlf%s.fasta' % (os.getcwd(), donor_gene, n)
        donor_r = '%s/%s_drf%s.fasta' % (os.getcwd(), donor_gene, n)
        output_1 = ''
        output_2 = ''
        command_blast_1 = ''
        command_blast_2 = ''
        if gene_r_strand == gene_d_strand:
            output_1 = '%s/%s_rlf%s___%s_dlf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            output_2 = '%s/%s_rrf%s___%s_drf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            command_blast_1 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_l, donor_l, output_1, blast_parameters)
            command_blast_2 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_r, donor_r, output_2, blast_parameters)
        if gene_r_strand != gene_d_strand:
            output_1 = '%s/%s_rlf%s___%s_drf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            output_2 = '%s/%s_rrf%s___%s_dlf%s.tab' % (os.getcwd(), recipient_gene, n, donor_gene, n)
            command_blast_1 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_l, donor_r, output_1, blast_parameters)
            command_blast_2 = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, query_r, donor_l, output_2, blast_parameters)
        os.system(command_blast_1)
        os.system(command_blast_2)

        # get recipient_left_iden_profile
        if os.path.isfile(output_1) and (os.path.getsize(output_1) == 0):
            recipient_left_iden_profile.append(['None'])
        if os.path.getsize(output_1) > 0:
            identity_l_list = []
            for each_line_1 in open(output_1):
                each_line_1_split = each_line_1.strip().split('\t')
                identity_1 = round(float(each_line_1_split[2]))
                alignment_length_1 = int(each_line_1_split[3])
                if alignment_length_1 >= alignment_length_cutoff:
                    identity_l_list.append(identity_1)
            if identity_l_list == []:
                recipient_left_iden_profile.append(['None'])
            else:
                recipient_left_iden_profile.append(identity_l_list)

        # get recipient_right_iden_profile
        if os.path.isfile(output_2) and (os.path.getsize(output_2) == 0):
            recipient_right_iden_profile.append(['None'])
        if os.path.getsize(output_2) > 0:
            identity_r_list = []
            for each_line_2 in open(output_2):
                each_line_2_split = each_line_2.strip().split('\t')
                identity_2 = round(float(each_line_2_split[2]))
                alignment_length_2 = int(each_line_2_split[3])
                if alignment_length_2 >= alignment_length_cutoff:
                    identity_r_list.append(identity_2)
            if identity_r_list == []:
                recipient_right_iden_profile.append(['None'])
            else:
                recipient_right_iden_profile.append(identity_r_list)
        n += 1

    print(recipient_left_iden_profile)
    print(recipient_right_iden_profile)

    # analyze flanking profile
    recipient_left_iden_profile_uniq = uniq_list(recipient_left_iden_profile)
    recipient_right_iden_profile_uniq = uniq_list(recipient_right_iden_profile)

    match_category = ''
    if (len(recipient_left_iden_profile_uniq) == 1) and (len(recipient_right_iden_profile_uniq) == 1) and (recipient_left_iden_profile_uniq[0][0] == 'None') and (recipient_right_iden_profile_uniq[0][0] == 'None'):
        match_category = 'Uniq matched'
    elif (recipient_left_iden_profile[-1][0] != 'None') or (recipient_right_iden_profile[-1][0] != 'None'):
        match_category = 'End matched'
    elif ((len(recipient_left_iden_profile_uniq) > 1) or (len(recipient_right_iden_profile_uniq) > 1)) and ((recipient_left_iden_profile[-1][0] == 'None') and (recipient_right_iden_profile[-1][0] == 'None')):
        match_category = 'Non_end multiple matched'

    return match_category


# define input files
flanking_length = 3000
calculation_step = 1000
pwd_blastn_exe = 'blastn'

#folder_name = 'AAM_00209___BNM_00032'
#folder_name = 'AAM_00308___BRT_00724'
#folder_name = 'ASN_02461___BHS_02663'
#folder_name = 'AAM_00754___BAD_01373'
#folder_name = 'AMAC_01101___BAD_02659'
#folder_name = 'BRT_03687___ARS_03341'
folder_name = 'AA_Refined_8_00597___CF_Refined_23_00314'


match_category = get_match_category(folder_name, flanking_length, calculation_step, pwd_blastn_exe)
print('%s: %s' % (folder_name, match_category))



