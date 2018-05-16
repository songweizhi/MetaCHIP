from Bio import SeqIO

candidates_file = open('/Users/weizhisong/Desktop/0_ready_for_submit(quiver_polished)/keep_list.txt')

# get candidates list:
candidate_list = []
for each in candidates_file:
    each = each.strip().split(' ')[0]
    candidate_list.append(each)
print(candidate_list)

for seq_record in SeqIO.parse('/Users/weizhisong/Desktop/000/combined_10.consensus.fasta', 'fasta'):
    print(seq_record.id)
    for candidate in candidate_list:
        if seq_record.id == candidate:
            SeqIO.write(seq_record, '/Users/weizhisong/Desktop/000/' + candidate + '.fasta', 'fasta')

