import os
import glob
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


os.chdir('/Users/songweizhi/Desktop')


input_gbk_folder = '/Users/songweizhi/Desktop/input_gbk_files'
prefix = 'NorthSea'
output_folder = '%s_MetaCHIP_outputs' % prefix

# os.mkdir(output_folder)

# prefix
# output_folder
# input_gbk_folder


faa_folder = '%s_faa_files' % prefix
ffn_folder = '%s_ffn_files' % prefix
pwd_faa_folder = '%s/%s' % (output_folder, faa_folder)
pwd_ffn_folder = '%s/%s' % (output_folder, ffn_folder)


# os.mkdir(pwd_faa_folder)
# os.mkdir(pwd_ffn_folder)

# get input gbk file list
input_gbk_re = '%s/*.gbk' % input_gbk_folder
input_gbk_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_gbk_re)]

for gbk_file in input_gbk_file_list:

    # prepare file name
    gbk_file_basename, gbk_file_extension = os.path.splitext(gbk_file)
    pwd_gbk_file = '%s/%s' % (input_gbk_folder, gbk_file)
    pwd_output_ffn_file = '%s/%s.ffn' % (pwd_ffn_folder, gbk_file_basename)
    pwd_output_faa_file = '%s/%s.faa' % (pwd_faa_folder, gbk_file_basename)

    output_ffn_handle = open(pwd_output_ffn_file, 'w')
    output_faa_handle = open(pwd_output_faa_file, 'w')
    for seq_record in SeqIO.parse(pwd_gbk_file, 'genbank'):
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
                export_dna_record(seq_nc, feature_id, feature_description, output_ffn_handle)
                export_aa_record(seq_aa, feature_id, feature_description, output_faa_handle)

    output_ffn_handle.close()
    output_faa_handle.close()




