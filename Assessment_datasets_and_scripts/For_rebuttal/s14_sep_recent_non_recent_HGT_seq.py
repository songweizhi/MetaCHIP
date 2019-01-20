from Bio import SeqIO


seq_file =                  '/Users/songweizhi/Desktop/nature_rebuttal/Total2094_g664_HGTs_PG_aa.fasta'
PG_validated_HGTs =         '/Users/songweizhi/Desktop/nature_rebuttal/Total2094_g664_HGTs_PG_validated.txt'
seq_file_recent_HGT =       '/Users/songweizhi/Desktop/nature_rebuttal/Total2094_g664_HGTs_PG_recent.faa'
seq_file_non_recent_HGT =   '/Users/songweizhi/Desktop/nature_rebuttal/Total2094_g664_HGTs_PG_non_recent.faa'

recent_HGTs = set()
non_recent_HGTs = set()
for each in open(PG_validated_HGTs):
    each_split = each.strip().split('\t')

    if not each.startswith('Gene_1'):

        gene_1 = each_split[0]
        gene_2 = each_split[1]
        identity = float(each_split[4])

        if identity >= 99:
            recent_HGTs.add(gene_1)
            recent_HGTs.add(gene_2)
        else:
            non_recent_HGTs.add(gene_1)
            non_recent_HGTs.add(gene_2)


seq_file_recent_HGT_handle = open(seq_file_recent_HGT, 'w')
seq_file_non_recent_HGT_handle = open(seq_file_non_recent_HGT, 'w')
for seq in SeqIO.parse(seq_file, 'fasta'):

    if seq.id in recent_HGTs:
        seq_file_recent_HGT_handle.write('>%s\n' % seq.id)
        seq_file_recent_HGT_handle.write('%s\n' % str(seq.seq))

    if seq.id in non_recent_HGTs:
        seq_file_non_recent_HGT_handle.write('>%s\n' % seq.id)
        seq_file_non_recent_HGT_handle.write('%s\n' % str(seq.seq))

seq_file_recent_HGT_handle.close()
seq_file_non_recent_HGT_handle.close()

