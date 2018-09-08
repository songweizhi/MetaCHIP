import os
import glob
import shutil
import argparse
import subprocess
from time import sleep
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


# to do
# include filter with blast coverage?


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def get_best_domain_hit_seq(pwd_domtblout_file, protein_Sequence_dict, pwd_output_folder):

    current_hmm_id = ''
    current_target_gene = ''
    current_hmm_pos1 = 0
    current_hmm_pos2 = 0
    current_hmm_score = 0
    # pwd_best_domain_hit_file_handle = open(pwd_best_domain_hit_file, 'w')
    for each_hit in open(pwd_domtblout_file):
        if not each_hit.startswith('#'):
            each_hit_split = each_hit.strip().split(' ')
            each_hit_split_no_space = []
            for each_element in each_hit_split:
                if each_element != '':
                    each_hit_split_no_space.append(each_element)
            target_gene = each_hit_split_no_space[0]
            hmm_id = each_hit_split_no_space[4]
            hmm_score = float(each_hit_split_no_space[13])
            hmm_pos1 = int(each_hit_split_no_space[17])
            hmm_pos2 = int(each_hit_split_no_space[18])

            if current_hmm_id == '':
                current_hmm_id = hmm_id
                current_target_gene = target_gene
                current_hmm_pos1 = hmm_pos1
                current_hmm_pos2 = hmm_pos2
                current_hmm_score = hmm_score
            elif current_hmm_id == hmm_id:
                if hmm_score > current_hmm_score:
                    current_hmm_id = hmm_id
                    current_target_gene = target_gene
                    current_hmm_pos1 = hmm_pos1
                    current_hmm_pos2 = hmm_pos2
                    current_hmm_score = hmm_score
            elif current_hmm_id != hmm_id:
                seqout_file_basename = '%s___%s' % (seqin_file_name, current_hmm_id)
                seqout_file_name = '%s.faa' % (seqout_file_basename)
                pwd_seqout_file_name = '%s/%s' % (pwd_output_folder, seqout_file_name)
                seqout_file_handle = open(pwd_seqout_file_name, 'w')
                export_dna_record(protein_Sequence_dict[current_target_gene][(current_hmm_pos1 - 1):current_hmm_pos2], seqout_file_basename, '', seqout_file_handle)
                seqout_file_handle.close()
                # pwd_best_domain_hit_file_handle.write('%s\t%s\t%s\t%s\t%s\n' % (current_hmm_id, current_target_gene, current_hmm_score, current_hmm_pos1, current_hmm_pos2))
                current_hmm_id = hmm_id
                current_target_gene = target_gene
                current_hmm_pos1 = hmm_pos1
                current_hmm_pos2 = hmm_pos2
                current_hmm_score = hmm_score
    # pwd_best_domain_hit_file_handle.write('%s\t%s\t%s\t%s\t%s\n' % (current_hmm_id, current_target_gene, current_hmm_score, current_hmm_pos1, current_hmm_pos2))
    # pwd_best_domain_hit_file_handle.close()


def get_assignment(list_in, confidence_level):

    # uniq list
    list_in_uniq = []
    for each_taxa in list_in:
        if each_taxa not in list_in_uniq:
            list_in_uniq.append(each_taxa)

    # get percentage for each uniq taxa and calculate assignment
    assignment = ''
    assignment_percentage = 0
    for each_uniq_taxa in list_in_uniq:
        each_uniq_taxa_count = list_in.count(each_uniq_taxa)
        each_uniq_taxa_percent = each_uniq_taxa_count * 100 / len(list_in)
        each_uniq_taxa_percent = float("{0:.2f}".format(each_uniq_taxa_percent))
        each_uniq_taxa_split = each_uniq_taxa.split('__')
        if (each_uniq_taxa_percent >= confidence_level) and (each_uniq_taxa_split[1] != ''):
            assignment = each_uniq_taxa
            assignment_percentage = each_uniq_taxa_percent

    return assignment, assignment_percentage


######################################################## inputs ########################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i',
                        required=True,
                        help='folder holds the annotation file of input genomes')

    parser.add_argument('-c',
                        required=False,
                        type=int,
                        default=85,
                        help='confidence level, default: 85')

    parser.add_argument('-tuning',
                        required=False,
                        action="store_true",
                        help='tuning mode')

    parser.add_argument('-bac120hmm',
                        required=False,
                        default='/Users/songweizhi/Softwares/get_taxon_GTDB/bac120.HMM',
                        help='path to bac120.HMM')

    parser.add_argument('-blastdb',
                        required=False,
                        default='/Users/songweizhi/Softwares/get_taxon_GTDB/bac120_msa_individual_genes_r83_no_dash',
                        help='path to folder holds blast database: bac120_msa_individual_genes_r83_no_dash')

    parser.add_argument('-taxondb',
                        required=False,
                        default='/Users/songweizhi/Softwares/get_taxon_GTDB/taxonomy_r83_March2018.tsv',
                        help='path to file: taxonomy_r83_March2018.tsv')

    parser.add_argument('-hmmsearch',
                        required=False,
                        default='/Users/songweizhi/Softwares/hmmer/hmmer-3.1b2-macosx-intel/binaries/hmmsearch',
                        help='path to hmmsearch')

    parser.add_argument('-blastp',
                        required=False,
                        default='blastp',
                        help='path to blastp')

    parser.add_argument('-ktImportText',
                        required=False,
                        default='ktImportText',
                        help='path to ktImportText')

    parser.add_argument('-keep_tmp',
                        action="store_true",
                        required=False,
                        help='keep temporary files')

    args = vars(parser.parse_args())
    genome_folder = args['i']
    confidence_level = args['c']
    pwd_hmm_db = args['bac120hmm']
    pwd_blastdb = args['blastdb']
    gtdb_taxonomy_db = args['taxondb']
    pwd_hmmsearch_exe = args['hmmsearch']
    pwd_blastp_exe = args['blastp']
    pwd_ktImportText_executable = args['ktImportText']
    tuning_mode = args['tuning']
    keep_tmp = args['keep_tmp']

    if genome_folder[-1] == '/':
        input_bin_folder_1 = genome_folder[:-1]

    #################################################### define file name ##################################################

    wd = os.getcwd()

    op_folder =                     '%s_outputs' % genome_folder.split('/')[-1]
    op_folder_temp =                'temp'
    op_folder_krona =               'Krona_plots'
    op_folder_krona_input_txt =     'Krona_plots_input_txt'
    output_finest =                 'taxon_assignment.txt'
    output_lineage =                'taxon_assignment_lineage.txt'

    pwd_op_folder =                 '%s/%s'    % (wd, op_folder)
    pwd_op_folder_temp =            '%s/%s'    % (pwd_op_folder, op_folder_temp)
    pwd_op_folder_krona =           '%s/%s'    % (pwd_op_folder, op_folder_krona)
    pwd_op_folder_krona_input_txt = '%s/%s'    % (pwd_op_folder, op_folder_krona_input_txt)
    pwd_output_finest =             '%s/%s'    % (pwd_op_folder, output_finest)
    pwd_output_lineage =            '%s/%s'    % (pwd_op_folder, output_lineage)

    # create temporary folder
    if tuning_mode != True:
        if os.path.isdir(pwd_op_folder):
            os.system('rm -r %s' % pwd_op_folder)
            os.mkdir(pwd_op_folder)
            os.mkdir(pwd_op_folder_temp)
            os.mkdir(pwd_op_folder_krona)
            os.mkdir(pwd_op_folder_krona_input_txt)
        else:
            os.mkdir(pwd_op_folder)
            os.mkdir(pwd_op_folder_temp)
            os.mkdir(pwd_op_folder_krona)
            os.mkdir(pwd_op_folder_krona_input_txt)

    # get input genome file list
    input_genome_re = '%s/*.%s' % (genome_folder, 'faa')
    input_genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_re)]
    input_genome_file_list = sorted(input_genome_file_list)

    print()

    if len(input_genome_file_list) == 0:
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' No input genome detected, program exited')
        exit()

    sleep(0.5)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' The number of input genome: %s' % (len(input_genome_file_list)))
    sleep(0.5)

    pwd_output_finest_handle = open(pwd_output_finest, 'w')
    pwd_output_finest_handle.write('Genome\tTaxon\tPercentage\n')
    pwd_output_finest_handle.close()

    pwd_output_lineage_handle = open(pwd_output_lineage, 'w')
    pwd_output_lineage_handle.close()

    n = 1
    for each_query_genome in input_genome_file_list:
        pwd_query_genome = '%s/%s' % (genome_folder, each_query_genome)
        seqin_file_name, faa_file_ext = os.path.splitext(each_query_genome)

        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Processing %s/%s: %s' % (n, len(input_genome_file_list), each_query_genome))

        ################################################## define file name ################################################

        output_subfolder =              '%s'                                % seqin_file_name
        domtblout_file =                '%s_hmmout.tsv'                     % seqin_file_name
        combined_blast_results =        '%s_combined_blast_results.txt'     % seqin_file_name
        identified_lineage =            '%s_identified_lineage.txt'         % seqin_file_name
        identified_lineage_sorted =     '%s_identified_lineage_sorted.txt'  % seqin_file_name
        krona_in =                      '%s_krona_in.txt'                   % seqin_file_name
        krona_out =                     '%s.html'                           % seqin_file_name
        krona_report =                  '%s_krona_report.txt'               % seqin_file_name

        pwd_output_subfolder =          '%s/%s'                     % (pwd_op_folder_temp, output_subfolder)
        pwd_domtblout_file =            '%s/%s'                     % (pwd_output_subfolder, domtblout_file)
        pwd_combined_blast_results =    '%s/%s'                     % (pwd_output_subfolder, combined_blast_results)
        pwd_identified_lineage =        '%s/%s'                     % (pwd_output_subfolder, identified_lineage)
        pwd_identified_lineage_sorted = '%s/%s'                     % (pwd_output_subfolder, identified_lineage_sorted)
        pwd_krona_in =                  '%s/%s'                     % (pwd_op_folder_krona_input_txt, krona_in)
        pwd_krona_out =                 '%s/%s'                     % (pwd_op_folder_krona, krona_out)
        pwd_krona_report =              '%s/%s'                     % (pwd_output_subfolder, krona_report)

        ####################################################################################################################

        # create temporary folder
        if tuning_mode != True:
            if os.path.isdir(pwd_output_subfolder):
                shutil.rmtree(pwd_output_subfolder)
                os.makedirs(pwd_output_subfolder)
            else:
                os.makedirs(pwd_output_subfolder)

        # get protein_Sequence_dict
        protein_Sequence_dict = {}
        for seq_record in SeqIO.parse(pwd_query_genome, 'fasta'):
            protein_Sequence_dict[seq_record.id] = str(seq_record.seq)

        # run hmmsearch
        if tuning_mode != True:
            hmmsearch_cmd = '%s -o /dev/null --domtblout %s %s %s' % (pwd_hmmsearch_exe, pwd_domtblout_file, pwd_hmm_db, pwd_query_genome)
            os.system(hmmsearch_cmd)

        # get best domain hit sequences
        get_best_domain_hit_seq(pwd_domtblout_file, protein_Sequence_dict, pwd_output_subfolder)

        # get file list
        file_re = '%s/*.%s' % (pwd_output_subfolder, 'faa')
        file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]
        file_list = sorted(file_list)
        MaxProcesses = 20
        Processes = []

        for each_file in file_list:
            each_file_name, each_file_ext = os.path.splitext(each_file)
            hmm_id = each_file_name.split('___')[1]
            pwd_each_file = '%s/%s' % (pwd_output_subfolder, each_file)
            pwd_each_file_db = '%s/gtdb_%s.faa' % (pwd_blastdb, hmm_id)
            pwd_each_file_blastp_out = '%s/%s_blastp.txt' % (pwd_output_subfolder, each_file_name)
            blastp_cmd = '%s -query %s -subject %s -outfmt 6 -out %s -max_target_seqs 1' % (pwd_blastp_exe, pwd_each_file, pwd_each_file_db, pwd_each_file_blastp_out)

            # run without subprocess
            # if tuning_mode != True:
            #     os.system(blastp_cmd)

            # run with subprocess
            blastp_cmd_list = [pwd_blastp_exe, '-query', pwd_each_file, '-subject', pwd_each_file_db, '-outfmt', '6', '-out', pwd_each_file_blastp_out, '-max_target_seqs', '1']
            if tuning_mode != True:
                # keep wait if there is no spare slots
                while len(Processes) >= MaxProcesses:
                    sleep(0.1)
                    for process in Processes:
                        if process.poll() is not None:
                            Processes.remove(process)
                # submit new subprocess
                p = subprocess.Popen(blastp_cmd_list)
                Processes.append(p)
        # wait for all subprocess to finish
        if tuning_mode != True:
            p.wait()

        # combine blast results together
        if tuning_mode != True:
            os.system('cat %s/*_blastp.txt > %s' % (pwd_output_subfolder, pwd_combined_blast_results))

        taxon_id_to_lineage_dict = {}
        for each_taxon in open(gtdb_taxonomy_db):
            each_taxon_split = each_taxon.strip().split('\t')
            each_taxon_id = each_taxon_split[0]
            each_taxon_lineage = each_taxon_split[1]
            taxon_id_to_lineage_dict[each_taxon_id] = each_taxon_lineage

        # get reference id list
        identified_lineage_handle = open(pwd_identified_lineage, 'w')
        for each in open(pwd_combined_blast_results):
            each_split = each.strip().split('\t')
            ref_id = each_split[1]
            identified_lineage_handle.write('%s\n' % (taxon_id_to_lineage_dict[ref_id]))
        identified_lineage_handle.close()

        ################################################## get Krona plot ##################################################

        # prepare input for Krona plot
        os.system('cat %s | sort > %s' % (pwd_identified_lineage, pwd_identified_lineage_sorted))

        identified_lineage_list = []
        identified_lineage_list_uniq = []
        for each_identified_lineage in open(pwd_identified_lineage_sorted):
            each_identified_lineage = each_identified_lineage.strip()
            identified_lineage_list.append(each_identified_lineage)
            if each_identified_lineage not in identified_lineage_list_uniq:
                identified_lineage_list_uniq.append(each_identified_lineage)

        identified_lineage_uniq_count_dict = {}
        for each_identified_lineage_uniq in identified_lineage_list_uniq:
            current_count = identified_lineage_list.count(each_identified_lineage_uniq)
            identified_lineage_uniq_count_dict[each_identified_lineage_uniq] = current_count

        # get input file for Krona plot
        pwd_krona_in_handle = open(pwd_krona_in, 'w')
        for each_key in identified_lineage_uniq_count_dict:
            each_key_split = each_key.split(';')
            pwd_krona_in_handle.write('%s\t%s\n' % (identified_lineage_uniq_count_dict[each_key], '\t'.join(each_key_split)))
        pwd_krona_in_handle.close()

        # plot Krona
        os.system('%s %s -o %s > %s' % (pwd_ktImportText_executable, pwd_krona_in, pwd_krona_out, pwd_krona_report))

        ####################################################################################################################

        # get the number of assignment at all seven levels
        d_list = []
        p_list = []
        c_list = []
        o_list = []
        f_list = []
        g_list = []
        s_list = []
        for each_lineage in open(pwd_identified_lineage):
            each_lineage_split = each_lineage.strip().split(';')
            d_list.append(each_lineage_split[0])
            p_list.append(each_lineage_split[1])
            c_list.append(each_lineage_split[2])
            o_list.append(each_lineage_split[3])
            f_list.append(each_lineage_split[4])
            g_list.append(each_lineage_split[5])
            s_list.append(each_lineage_split[6])

        # put all list together for interation
        list_of_rank_list = [d_list, p_list, c_list, o_list, f_list, g_list, s_list]

        # get assignments at all rank levels
        assignment_list = []
        assignment_percentage_list = []
        assignment_combined = []
        for each_rank_list in list_of_rank_list:
            current_assignment, current_assignment_percentage = get_assignment(each_rank_list, confidence_level)
            if current_assignment_percentage >= confidence_level:
                assignment_list.append(current_assignment)
                assignment_percentage_list.append(current_assignment_percentage)
                combined = '%s(%s)' % (current_assignment, current_assignment_percentage)
                assignment_combined.append(combined)

        # get finest assignment
        finest_assignment = assignment_list[-1]
        finest_assignment_percent = assignment_percentage_list[-1]

        # write out
        pwd_output_finest_handle = open(pwd_output_finest, 'a')
        pwd_output_finest_handle.write('%s\t%s\t%s\n' % (seqin_file_name, finest_assignment, finest_assignment_percent))
        pwd_output_finest_handle.close()
        pwd_output_lineage_handle = open(pwd_output_lineage, 'a')
        pwd_output_lineage_handle.write('%s\t%s\n' % (seqin_file_name, ';'.join(assignment_combined)))
        pwd_output_lineage_handle.close()

        n += 1

    if keep_tmp != True:
        if os.path.isdir(pwd_op_folder_temp):
            shutil.rmtree(pwd_op_folder_temp, ignore_errors=True)
            if os.path.isdir(pwd_op_folder_temp):
                shutil.rmtree(pwd_op_folder_temp, ignore_errors=True)

    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' All done!')

