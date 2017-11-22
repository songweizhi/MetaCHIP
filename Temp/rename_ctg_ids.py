import os
import re
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


wd = '/Users/songweizhi/Desktop/rename_ctg_id'
input_bins = '%s/43bins_original/*.fasta' % (wd)
input_bin_list = [os.path.basename(file_name) for file_name in glob.glob(input_bins)]


n = 1
pwd_rename_log_handle = open('%s/rename_log.txt' % wd, 'w')
for each_bin in input_bin_list:
    print('Processing %s' % each_bin)
    each_bin_no_extension, ext = os.path.splitext(each_bin)
    pwd_each_bin = '%s/43bins_original/%s' % (wd, each_bin)
    pwd_each_bin_renamed = '%s/43bins_renamed/%s' % (wd, each_bin)
    pwd_each_bin_renamed_handle = open(pwd_each_bin_renamed, 'w')
    for each_ctg in SeqIO.parse(pwd_each_bin, 'fasta'):
        each_ctg_seq = str(each_ctg.seq)
        each_ctg_id = each_ctg.id
        new_id = 'scaffold_%s' % n
        export_dna_record(each_ctg_seq, new_id, '', pwd_each_bin_renamed_handle)
        pwd_rename_log_handle.write('%s\t%s\t%s\n' % (each_bin_no_extension, new_id, each_ctg_id))
        n += 1
    pwd_each_bin_renamed_handle.close()
pwd_rename_log_handle.close()
