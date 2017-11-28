
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


output_handle = open('/Users/songweizhi/Desktop/combined_16s.fasta', 'w')
for each_seq in SeqIO.parse('/Users/songweizhi/Desktop/combined.fasta', 'fasta'):
    seq_id = each_seq.id
    seq_id_new = '%s_16s_rRNA' % seq_id.split(':')[0]
    seq_sequence = str((each_seq.seq))
    if 2000 > len(each_seq.seq) > 1000:
        export_dna_record(seq_sequence, seq_id_new, '', output_handle)
output_handle.close()
