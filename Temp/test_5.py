import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def gbk2fasta(pwd_gbk_file):

    wd, gbk_file = os.path.split(pwd_gbk_file)
    gbk_file_name, gbk_file_ext = os.path.splitext(gbk_file)
    pwd_fasta_file = '%s/%s.fasta' % (wd, gbk_file_name)
    fasta_file = open(pwd_fasta_file, 'w')
    input_gbk = SeqIO.parse(pwd_gbk_file,'genbank')

    for record in input_gbk:
        new_record = SeqRecord(record.seq, id=record.id, description='')
        SeqIO.write(new_record, fasta_file, 'fasta')

gbk2fasta('/Users/songweizhi/Desktop/test.gbk')
