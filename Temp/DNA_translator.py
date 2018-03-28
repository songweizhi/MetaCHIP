from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna


output_protein_handle = open('/Users/songweizhi/Desktop/HGT_candidates_aa.fasta', 'w')
for each in SeqIO.parse('/Users/songweizhi/Desktop/HGT_candidates.fasta', 'fasta'):
    print(each.seq)
    aa = each.seq.translate()
    print(aa)
    print('\n')
    output_protein_handle.write('>%s %s\n%s\n' % (each.id, '', aa))


output_protein_handle.close()

