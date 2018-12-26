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


files = '/Users/songweizhi/Desktop/BetterBins/*.fa'
file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]
print(file_list)
print(len(file_list))



for each_bin in file_list:
    bin_file_basename = '.'.join(each_bin.split('.')[:-1])
    pwd_bin = '/Users/songweizhi/Desktop/BetterBins/%s' % each_bin
    pwd_bin_good = '/Users/songweizhi/Desktop/BetterBins_good/%s' % each_bin
    pwd_bin_good_handle = open(pwd_bin_good, 'w')

    bin_total_len = 0
    n = 1
    for ctg in SeqIO.parse(pwd_bin, 'fasta'):
        ctg_seq = str(ctg.seq)
        if 'N' in ctg_seq:
            ctg_seq_split = ctg_seq.split('N')
            for ctg_frag in ctg_seq_split:
                if len(ctg_frag) >= 2000:
                    ctg_frag_id = '%s_ctg%s' % (bin_file_basename, n)
                    export_dna_record(ctg_frag, ctg_frag_id, '', pwd_bin_good_handle)
                    bin_total_len += len(ctg_frag)
                    n += 1
        else:
            if len(ctg_seq) >= 2000:
                ctg_id_new = '%s_ctg%s' % (bin_file_basename, n)
                export_dna_record(ctg_seq, ctg_id_new, '', pwd_bin_good_handle)
                bin_total_len += len(ctg_seq)
                n += 1

    pwd_bin_good_handle.close()

    if bin_total_len < 204800:
        os.system('rm %s' % pwd_bin_good)





