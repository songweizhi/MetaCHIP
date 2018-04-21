import os
import glob
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from datetime import datetime

#
# def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
#     seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
#     seq_record = SeqRecord(seq_object)
#     seq_record.id = gene_id
#     seq_record.description = gene_description
#     SeqIO.write(seq_record, output_handle, 'fasta')
#
#
# # get file list
# file_re = '/Users/songweizhi/Desktop/test_hmm/bac120_msa_individual_genes_r83_no_dash/*.faa'
# file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]
# file_list = sorted(file_list)
#
# combined_faa_db = open('/Users/songweizhi/Desktop/test_hmm/combined_db.faa', 'w')
# for each_file in file_list:
#     pwd_eachfile = '/Users/songweizhi/Desktop/test_hmm/bac120_msa_individual_genes_r83_no_dash/%s' % each_file
#     each_file_name, each_file_ext = os.path.splitext(each_file)
#     for each_seq in SeqIO.parse(pwd_eachfile, 'fasta'):
#         new_id = '%s___%s' % (each_file_name, each_seq.id)
#         export_dna_record(str(each_seq.seq), new_id, '', combined_faa_db)
#
# combined_faa_db.close()


os.chdir('/Users/songweizhi/Desktop/test_hmm')

print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
os.system('blastp -query combined_bin1216.faa -subject combined_db.faa -outfmt 6 -out combined_bin1216_blast.txt -max_target_seqs 1')
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
