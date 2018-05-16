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


os.chdir('/Users/songweizhi/Desktop/test_end_break')
input_file = 'combined_ctgs.fasta'
end_length = 5000
output_file = 'combined_ctgs_%sbp.fasta' % end_length

output_handle = open(output_file, 'w')
for each_seq_record in SeqIO.parse(input_file, 'fasta'):
    sequence = str(each_seq_record.seq)
    sequence_left = sequence[0:end_length]
    sequence_left_id = '%s_l%sbp' % (each_seq_record.id, end_length)
    sequence_right = sequence[-end_length:]
    sequence_right_id = '%s_r%sbp' % (each_seq_record.id, end_length)
    if len(sequence) > 2 * end_length:
        export_dna_record(sequence_left, sequence_left_id, '', output_handle)
        export_dna_record(sequence_right, sequence_right_id, '', output_handle)
output_handle.close()


# Run Blast
path_to_makeblastdb_executable = 'makeblastdb'
makeblastdb_cmd = path_to_makeblastdb_executable + ' -in %s -dbtype nucl -parse_seqids' % output_file
os.system(makeblastdb_cmd)
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'

blast_cmd = 'blastn -query %s -db %s -out all_vs_all_%sbp.tab %s' % (output_file, output_file, end_length, blast_parameters)

os.system(blast_cmd)


