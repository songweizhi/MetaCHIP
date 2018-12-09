import os
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


wd = '/Users/songweizhi/Desktop/GGG2'
query_seq_file = 'gene_seq_bootstrap_2.fasta'
subject_seq_file = 'bootstrap2_m0_combined.ffn'
blast_output = 'blast_out2.tab'



os.chdir(wd)
#
# blastn_cmd = 'blastn -query %s -subject %s -outfmt 6 -out %s -task blastn' % (query_seq_file, subject_seq_file, blast_output)
# os.system(blastn_cmd)

# blastn -query 53054_purH.fna -subject 53054_purH.fna -outfmt 6 -out 53054_purH.tab
# blastn -query 53054_purH.fna -subject 53054_purH.fna -outfmt 6 -out 53054_purH.tab -task blastn
# blastn -query AKV_01260.fasta -subject bootstrap2_m0_combined.ffn -outfmt 6 -out AKV_01260.tab -task blastn
# blastn -query AKV_01260.fasta -subject bootstrap2_m0_combined.ffn -outfmt 6 -out AKV_01260.tab -task blastn




for each in open('AKV_01260.tab'):
    #print(each)
    if int(each.split('\t')[3]) >= 200:
        print(each.strip())
