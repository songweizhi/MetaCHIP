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


# gz_file_folder = '/Users/songweizhi/Desktop/'
# renamed_file_folder = '/Users/songweizhi/Desktop/renamed_genomes'
# taxon_id_to_fna_file = '/Users/songweizhi/Desktop/nature_rebuttal/taxon_id_to_fna.txt'
gz_file_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes'
renamed_file_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes_renamed'
taxon_id_to_fna_file = '/srv/scratch/z5039045/MetaCHIP_rebuttal/taxon_id_to_fna.txt'


taxon_id_to_fna_file_dict = {}
gz_id_list = []
for each in open(taxon_id_to_fna_file):
    each_split = each.strip().split('\t')
    taxon_id_to_fna_file_dict[each_split[1]] = each_split[0]
    gz_id_list.append(each_split[1])


for file_id in gz_id_list:
    pwd_seq_file = '%s/%s_genomic.fna' % (gz_file_folder, file_id)
    new_file_name = 'taxon%s.fasta' % taxon_id_to_fna_file_dict[file_id]
    pwd_new_file = '%s/%s' % (renamed_file_folder, new_file_name)
    new_file_handle = open(pwd_new_file, 'w')
    ctg_index = 1
    for ctg in SeqIO.parse(pwd_seq_file, 'fasta'):
        ctg_id_new = 'taxon%sCtg%s' % (taxon_id_to_fna_file_dict[file_id], ctg_index)
        export_dna_record(str(ctg.seq), ctg_id_new, '', new_file_handle)
        ctg_index += 1
    new_file_handle.close()
