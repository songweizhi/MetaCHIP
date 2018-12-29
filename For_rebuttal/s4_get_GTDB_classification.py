from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


# gz_file_folder = '/Users/songweizhi/Desktop/'
# renamed_file_folder = '/Users/songweizhi/Desktop/renamed_genomes'
# taxon_id_to_fna_file = '/Users/songweizhi/Desktop/nature_rebuttal/taxon_id_to_fna.txt'
# gz_file_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes'
# renamed_file_folder = '/srv/scratch/z5039045/MetaCHIP_rebuttal/downloaded_genomes_renamed'
taxon_id_to_fna_file = '/Users/songweizhi/Desktop/nature_rebuttal/taxon_id_to_fna.txt'


taxon_to_ncbi_dict = {}
ncbi_id_list = []
for each in open(taxon_id_to_fna_file):
    each_split = each.strip().split('\t')
    taxon_id = each_split[0]
    ncbi_id = each_split[1][:15]
    taxon_to_ncbi_dict[taxon_id] = ncbi_id
    ncbi_id_list.append(ncbi_id)

print(taxon_to_ncbi_dict)
print(ncbi_id_list)


taxonomy_r86 = '/Users/songweizhi/Desktop/nature_rebuttal/taxonomy_r86_August2018.tsv'


m = 0
for each in open(taxonomy_r86):
    ncbi_id = each.strip().split('\t')[0][3:]
    if ncbi_id in ncbi_id_list:
        print(each.strip())
        m += 1

print(m)
