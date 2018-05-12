import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def export_aa_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.protein)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


input_gbk_file = '/Users/songweizhi/Desktop/aaa/NorthSea_bin019.gbk'
export_type = 'A'  # N or P or A


# prepare file name
gbk_file_path, gbk_file_name = os.path.split(input_gbk_file)
gbk_file_name_no_extension, gbk_file_extension = os.path.splitext(gbk_file_name)
output_ffn_filename = '%s.ffn' % gbk_file_name_no_extension
output_faa_filename = '%s.faa' % gbk_file_name_no_extension

if gbk_file_path != '':
    pwd_output_ffn_file = '%s/%s' % (gbk_file_path, output_ffn_filename)
    pwd_output_faa_file = '%s/%s' % (gbk_file_path, output_faa_filename)
else:
    pwd_output_ffn_file = output_ffn_filename
    pwd_output_faa_file = output_faa_filename

output_ffn_handle = None
output_faa_handle = None
if export_type == 'N':
    output_ffn_handle = open(pwd_output_ffn_file, 'w')
if export_type == 'P':
    output_faa_handle = open(pwd_output_faa_file, 'w')
if export_type == 'A':
    output_ffn_handle = open(pwd_output_ffn_file, 'w')
    output_faa_handle = open(pwd_output_faa_file, 'w')

for seq_record in SeqIO.parse(input_gbk_file, 'genbank'):
    for feature in seq_record.features:
        if feature.type == "CDS":
            seq_record_sequence = str(seq_record.seq)

            # get DNA sequence
            seq_nc = ''
            if feature.location.strand == 1:
                seq_nc = seq_record_sequence[feature.location.start:feature.location.end]
            if feature.location.strand == -1:
                nc_seq_rc = seq_record_sequence[feature.location.start:feature.location.end]
                seq_nc = str(Seq(nc_seq_rc, generic_dna).reverse_complement())

            # get aa sequence
            seq_aa = feature.qualifiers['translation'][0]
            feature_id = feature.qualifiers['locus_tag'][0]
            feature_description = feature.qualifiers['product'][0]

            # export to file
            if export_type == 'N':
                export_dna_record(seq_nc, feature_id, feature_description, output_ffn_handle)
            if export_type == 'P':
                export_aa_record(seq_aa, feature_id, feature_description, output_faa_handle)
            if export_type == 'A':
                export_dna_record(seq_nc, feature_id, feature_description, output_ffn_handle)
                export_aa_record(seq_aa, feature_id, feature_description, output_faa_handle)

# close file
if export_type == 'N':
    output_ffn_handle.close()
if export_type == 'P':
    output_faa_handle.close()
if export_type == 'A':
    output_ffn_handle.close()
    output_faa_handle.close()




