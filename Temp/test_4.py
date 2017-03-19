from Bio import SeqIO


path_to_gbk_file = '/Users/songweizhi/Desktop/AAM.gbk'

records = SeqIO.parse(path_to_gbk_file, 'genbank')

for record in records:
    print(record.features)
    print(type(record.features))
    print(len(record.features))
    all_records = record.features
    all_records_new = []
    for each_gene in record.features:
        if 'locus_tag' in each_gene.qualifiers:

            if each_gene.qualifiers['locus_tag'][0] == 'AAM_03457':
                record.features.pop(each_gene)

                print(each_gene)
    print(record.features)
    print(len(record.features))




