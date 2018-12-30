import os
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


# get file list
files = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05/*.fasta'
#files = '/Users/songweizhi/Desktop/111/*.fasta'

file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]
print(file_list)

out_folder = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05_renamed'

for each_file in file_list:

    pwd_each_file = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05/%s' % each_file
    pwd_each_file_out = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05_renamed/%s' % each_file

    pwd_each_file_out_handle = open(pwd_each_file_out, 'w')
    for each_seq in SeqIO.parse(pwd_each_file, 'fasta'):
        new_id = '_'.join(each_seq.id.split('_')[:2])
        new_description = '_'.join(each_seq.id.split('_')[2:])
        export_dna_record(str(each_seq.seq), new_id, new_description, pwd_each_file_out_handle)

    pwd_each_file_out_handle.close()




