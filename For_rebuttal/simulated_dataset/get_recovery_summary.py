import os
from Bio import SeqIO

usage = '''

module load python/3.5.2
module load blast+/2.6.0
python3 get_recovery_summary.py

'''


def get_recovered_HGTs_PC(gene_distribution, transferred_gene_seq, donor_genome_marker, recipient_genome_marker, HGT_BM_seq, HGT_PG_seq, HGT_candidates_PG_validated):

    blastn_BM_result = 'blastn_BM.tab'
    blastn_PG_result = 'blastn_PG.tab'
    blastn_BM = 'blastn -query %s -subject %s -outfmt 6 -out %s' % (transferred_gene_seq, HGT_BM_seq, blastn_BM_result)
    blastn_PG = 'blastn -query %s -subject %s -outfmt 6 -out %s' % (transferred_gene_seq, HGT_PG_seq, blastn_PG_result)
    os.system(blastn_BM)
    os.system(blastn_PG)

    # get the length of transferred genes
    gene_len_dict = {}
    for gene_seq in SeqIO.parse(transferred_gene_seq, 'fasta'):
        gene_len_dict[str(gene_seq.id)] = len(gene_seq.seq)

    # get gene flow direction dict
    gene_flow_direction_dict = {}
    for recipient in open(gene_distribution):
        recipient_split = recipient.strip().split(',')
        recipient_genome = recipient_split[0]
        transferred_genes = recipient_split[1:]
        for transferred_gene in transferred_genes:
            transferred_gene_genome = '_'.join(transferred_gene.split('_')[:-1])
            gene_flow_direction = '%s-->%s' % (transferred_gene_genome, recipient_genome)
            gene_flow_direction_dict[transferred_gene] = gene_flow_direction

    # get the number of BM recovered HGTs
    query_with_identical_HGT_BM = set()
    # print(blastn_BM_result)
    # os.system('cat %s' % blastn_BM_result)
    for BM_hit in open(blastn_BM_result):
        BM_hit_split = BM_hit.strip().split('\t')
        query_gene = BM_hit_split[0]
        subject_gene = BM_hit_split[1]
        identity = float(BM_hit_split[2])
        aln_len = int(BM_hit_split[3])
        if (subject_gene.startswith(recipient_genome_marker)) and (identity == 100) and ((aln_len / gene_len_dict[query_gene]) >= 0.98):
            query_with_identical_HGT_BM.add(query_gene)
    #print(query_with_identical_HGT_BM)
    # get the number of PG validated HGTs
    query_with_identical_HGT_PG = set()
    HGT_id_to_introduced_ID_dict = {}
    for PG_hit in open(blastn_PG_result):
        PG_hit_split = PG_hit.strip().split('\t')
        query_gene = PG_hit_split[0]
        subject_gene = PG_hit_split[1]
        identity = float(PG_hit_split[2])
        aln_len = int(PG_hit_split[3])
        if (subject_gene.startswith(recipient_genome_marker)) and (identity == 100) and ((aln_len / gene_len_dict[query_gene]) >= 0.98):
            query_with_identical_HGT_PG.add(query_gene)
            if subject_gene.startswith(recipient_genome_marker):
                HGT_id_to_introduced_ID_dict[subject_gene] = query_gene

    # get the number of PG validated HGTs with right direction
    HGT_with_right_direction = set()
    for PG_validated in open(HGT_candidates_PG_validated):
        PG_validated_split = PG_validated.strip().split('\t')
        gene_1 = PG_validated_split[0]
        gene_2 = PG_validated_split[1]
        gene_flow = PG_validated_split[7]

        gene_id_in_recipient = ''
        if gene_1.startswith(recipient_genome_marker):
            gene_id_in_recipient = gene_1
        else:
            gene_id_in_recipient = gene_2

        if (gene_id_in_recipient in HGT_id_to_introduced_ID_dict) and ((gene_flow_direction_dict[HGT_id_to_introduced_ID_dict[gene_id_in_recipient]]) == gene_flow):
            HGT_with_right_direction.add(gene_id_in_recipient)

    os.system('rm %s' % blastn_BM_result)
    os.system('rm %s' % blastn_PG_result)

    return len(query_with_identical_HGT_BM), len(query_with_identical_HGT_PG), len(HGT_with_right_direction)


def get_recovered_HGTs_PC_m0(gene_distribution, transferred_gene_seq, donor_genome_marker, recipient_genome_marker, HGT_BM_seq, HGT_PG_seq, HGT_candidates_PG_validated):

    blastn_BM_result = 'blastn_BM.tab'
    blastn_PG_result = 'blastn_PG.tab'
    blastn_BM = 'blastn -query %s -subject %s -outfmt 6 -out %s' % (transferred_gene_seq, HGT_BM_seq, blastn_BM_result)
    blastn_PG = 'blastn -query %s -subject %s -outfmt 6 -out %s' % (transferred_gene_seq, HGT_PG_seq, blastn_PG_result)
    os.system(blastn_BM)
    os.system(blastn_PG)

    # get the length of transferred genes
    gene_len_dict = {}
    for gene_seq in SeqIO.parse(transferred_gene_seq, 'fasta'):
        gene_len_dict[str(gene_seq.id)] = len(gene_seq.seq)

    # get gene flow direction dict
    gene_flow_direction_dict = {}
    for recipient in open(gene_distribution):
        recipient_split = recipient.strip().split(',')
        recipient_genome = recipient_split[0]
        transferred_genes = recipient_split[1:]
        for transferred_gene in transferred_genes:
            transferred_gene_genome = '_'.join(transferred_gene.split('_')[:-1])
            gene_flow_direction = '%s-->%s' % (transferred_gene_genome, recipient_genome)
            gene_flow_direction_dict[transferred_gene] = gene_flow_direction

    # get the number of BM recovered HGTs
    query_with_identical_HGT_BM = set()
    for BM_hit in open(blastn_BM_result):
        BM_hit_split = BM_hit.strip().split('\t')
        query_gene = BM_hit_split[0]
        subject_gene = BM_hit_split[1]
        identity = float(BM_hit_split[2])
        aln_len = int(BM_hit_split[3])
        if (identity == 100) and ((aln_len / gene_len_dict[query_gene]) >= 0.98):
            query_with_identical_HGT_BM.add(query_gene)

    # get the number of PG validated HGTs
    query_with_identical_HGT_PG = set()
    HGT_id_to_introduced_ID_dict = {}
    for PG_hit in open(blastn_PG_result):
        PG_hit_split = PG_hit.strip().split('\t')
        query_gene = PG_hit_split[0]
        subject_gene = PG_hit_split[1]
        identity = float(PG_hit_split[2])
        aln_len = int(PG_hit_split[3])
        if (identity == 100) and ((aln_len / gene_len_dict[query_gene]) >= 0.98):
            query_with_identical_HGT_PG.add(query_gene)
            if subject_gene.startswith(donor_genome_marker):
                HGT_id_to_introduced_ID_dict[subject_gene] = query_gene

    # get the number of PG validated HGTs with right direction
    HGT_with_right_direction = set()
    for PG_validated in open(HGT_candidates_PG_validated):
        PG_validated_split = PG_validated.strip().split('\t')
        gene_1 = PG_validated_split[0]
        gene_2 = PG_validated_split[1]
        gene_flow = PG_validated_split[7]

        gene_id_in_donor = ''
        if gene_1.startswith(donor_genome_marker):
            gene_id_in_donor = gene_1
        else:
            gene_id_in_donor = gene_2

        if (gene_id_in_donor in HGT_id_to_introduced_ID_dict) and (
                gene_flow_direction_dict[HGT_id_to_introduced_ID_dict[gene_id_in_donor]] == gene_flow):
            HGT_with_right_direction.add(HGT_id_to_introduced_ID_dict[gene_id_in_donor])

    os.system('rm %s' % blastn_BM_result)
    os.system('rm %s' % blastn_PG_result)

    return len(query_with_identical_HGT_BM), len(query_with_identical_HGT_PG), len(HGT_with_right_direction)


# cd to between group
donor_genome_marker = 'SB'
recipient_genome_marker = 'SM'
#donor_genome_marker = 'A'
#recipient_genome_marker = 'G'
mutation_levels = [0, 5, 10, 15, 20, 25, 30]


bootstrap_num = 1
while bootstrap_num <= 10:
    for mutation_level in mutation_levels:
        transferred_gene_seq = 'HgtSIM_wd/bootstrap%s/bootstrap%s_outputs_%s_1-0-1-1/input_sequence_mutant_nc.fasta' % (bootstrap_num, bootstrap_num, mutation_level)
        gene_distribution = 'HgtSIM_wd/HgtSIM_input_files/gene_distribution_bootstrap_%s.txt' % bootstrap_num


        # between class
        #HGT_BM_seq = 'MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_HGTs_ip90_al200bp_c70_e1000bp_f10kbp_c2/HGT_candidates_BM_nc.fasta' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        #HGT_PG_seq = 'MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_HGTs_ip90_al200bp_c70_e1000bp_f10kbp_c2/HGT_candidates_PG_nc.fasta' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        #HGT_candidates_PG_validated = 'MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_HGTs_ip90_al200bp_c70_e1000bp_f10kbp_c2/HGT_candidates_PG_validated.txt' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)


        # between genus
        HGT_BM_seq = 'MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_HGTs_ip90_al200bp_c70_e1000bp_f10kbp_g2/HGT_candidates_BM_nc.fasta' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        HGT_PG_seq = 'MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_HGTs_ip90_al200bp_c70_e1000bp_f10kbp_g2/HGT_candidates_PG_nc.fasta' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)
        HGT_candidates_PG_validated = 'MetaCHIP_wd/bootstrap%s/bootstrap%s_m%s/bootstrap%s_m%s_MetaCHIP_wd/bootstrap%s_m%s_HGTs_ip90_al200bp_c70_e1000bp_f10kbp_g2/HGT_candidates_PG_validated.txt' % (bootstrap_num, bootstrap_num, mutation_level, bootstrap_num, mutation_level, bootstrap_num, mutation_level)


        BM_recovered_HGTs = 0
        PG_recovered_HGTs = 0
        HGT_right_direction = 0
        if mutation_level == 0:
            BM_recovered_HGTs, PG_recovered_HGTs, HGT_right_direction = get_recovered_HGTs_PC_m0(gene_distribution, transferred_gene_seq, donor_genome_marker, recipient_genome_marker, HGT_BM_seq, HGT_PG_seq, HGT_candidates_PG_validated)
        else:
            BM_recovered_HGTs, PG_recovered_HGTs, HGT_right_direction = get_recovered_HGTs_PC(gene_distribution, transferred_gene_seq, donor_genome_marker, recipient_genome_marker, HGT_BM_seq, HGT_PG_seq, HGT_candidates_PG_validated)

        print('bootstrap%s_m%s\t%s\t%s\t%s' % (bootstrap_num, mutation_level, BM_recovered_HGTs, PG_recovered_HGTs, HGT_right_direction))

    bootstrap_num += 1















# wd = '/Users/songweizhi/Desktop/bootstrap1_m0'
# mutation_level = 0
# os.chdir(wd)
# transferred_gene_seq = 'input_sequence_mutant_nc.fasta'
# HGT_BM_seq = 'HGT_candidates_BM_nc.fasta'
# HGT_PG_seq = 'HGT_candidates_PG_nc.fasta'
# gene_distribution = 'gene_distribution_bootstrap_1.txt'
# HGT_candidates_PG_validated = 'HGT_candidates_PG_validated.txt'
#
#
#
# BM_recovered_HGTs = 0
# PG_recovered_HGTs = 0
# HGT_right_direction = 0
# if mutation_level == 0:
#     BM_recovered_HGTs, PG_recovered_HGTs, HGT_right_direction = get_recovered_HGTs_PC_m0(gene_distribution, transferred_gene_seq, donor_genome_marker, recipient_genome_marker, HGT_BM_seq, HGT_PG_seq, HGT_candidates_PG_validated)
# else:
#     BM_recovered_HGTs, PG_recovered_HGTs, HGT_right_direction = get_recovered_HGTs_PC(gene_distribution, transferred_gene_seq, donor_genome_marker, recipient_genome_marker, HGT_BM_seq, HGT_PG_seq, HGT_candidates_PG_validated)
#
#
# print('%s\t%s\t%s' % (BM_recovered_HGTs, PG_recovered_HGTs, HGT_right_direction))

