from Bio import SeqIO

path_to_gbk_file = '/Users/songweizhi/Desktop/AAM.gbk'
candidates = ['AAM_00001', 'AAM_00002', 'AAM_03075', 'AAM_03456', 'AAM_03457']

records = SeqIO.parse(path_to_gbk_file, 'genbank')
candidates_flanking_ranges = []
candidates_list = []
candidates_contig_list = []
candidates_start_list = []
candidates_end_list = []

for record in records:
    # get contig length
    contig_length = 0
    for each_gene in record.features:
        if each_gene.type == 'source':
            contig_length = int(each_gene.location.end)

        elif 'locus_tag' in each_gene.qualifiers:
            if each_gene.qualifiers["locus_tag"][0] in candidates:
               # get new start position
               each_gene_start = each_gene.location.start
               each_gene_start_50000 = each_gene_start - 50000
               if each_gene_start_50000 < 0:
                   each_gene_start_50000 = 0
               # get new end position
               each_gene_end = each_gene.location.end
               each_gene_end_50000 = each_gene.location.end + 50000
               if each_gene_end_50000 > contig_length:
                   each_gene_end_50000 = contig_length
               # append to list
               candidates_list.append(each_gene.qualifiers["locus_tag"][0])
               candidates_contig_list.append(record.id)
               candidates_start_list.append(each_gene_start_50000)
               candidates_end_list.append(each_gene_end_50000)
               candidates_flanking_ranges.append([each_gene.qualifiers["locus_tag"][0], record.id, each_gene_start_50000, each_gene_end_50000])


for each_candidate in candidates_flanking_ranges:
    can_id = each_candidate[0]
    can_contig_id = each_candidate[1]
    can_start = each_candidate[2]
    can_end = each_candidate[3]


records = SeqIO.parse(path_to_gbk_file, 'genbank')
for record in records:
    print(record)

