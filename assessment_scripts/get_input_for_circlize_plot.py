import os
import argparse

usage = """

python /Users/songweizhi/PycharmProjects/MetaCHIP/assessment_scripts/get_input_for_circlize_plot.py -in CF_HGT_candidates_ET_validated.txt 

"""

parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='input file')

args = vars(parser.parse_args())
input_file = args['in']

name2id_dict = {}
transfers = []
for each in open(input_file):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        Gene_1 = each_split[0]
        Gene_2 = each_split[1]
        Genome_1_ID = each_split[2]
        Genome_1 = '_'.join(Gene_1.split('_')[:-1])
        Genome_2_ID = each_split[3]
        Genome_2 = '_'.join(Gene_2.split('_')[:-1])
        Identity = each_split[4]
        End_break = each_split[5]
        Direction = each_split[6]
        if Genome_1 not in name2id_dict:
            name2id_dict[Genome_1] = Genome_1_ID
        if Genome_2 not in name2id_dict:
            name2id_dict[Genome_2] = Genome_2_ID
        transfers.append(Direction)

tmp1 = open('tmp1.txt', 'w')
all_group_id = []
for each_t in transfers:
    each_t_split = each_t.split('-->')
    donor = each_t_split[0]
    recipient = each_t_split[1]
    donor_id = name2id_dict[donor].split('_')[0]
    recipient_id = name2id_dict[recipient].split('_')[0]
    if donor_id not in all_group_id:
        all_group_id.append(donor_id)
    if recipient_id not in all_group_id:
        all_group_id.append(recipient_id)
    tmp1.write('%s,%s\n' % (donor_id, recipient_id))
tmp1.close()

os.system('cat tmp1.txt | sort > tmp1_sorted.txt')

current_t = ''
count = 0
tmp2 = open('tmp1_sorted_count.txt', 'w')
for each_t2 in open('tmp1_sorted.txt'):
    each_t2 = each_t2.strip()
    if current_t == '':
        current_t = each_t2
        count += 1
    elif current_t == each_t2:
        count += 1
    elif current_t != each_t2:
        tmp2.write('%s,%s\n' % (current_t, count))
        current_t = each_t2
        count = 1
tmp2.write('%s,%s\n' % (current_t, count))
tmp2.close()

# read in count as dict
transfer_count = {}
for each_3 in open('tmp1_sorted_count.txt'):
    each_3_split = each_3.strip().split(',')
    key = '%s,%s' % (each_3_split[0], each_3_split[1])
    value = each_3_split[2]
    transfer_count[key] = value

all_group_id = sorted(all_group_id)

input_file_name, input_file_ext = os.path.splitext(input_file)
matrix_file = open('%s_matrix.csv' % input_file_name, 'w')
matrix_file.write('\t' + '\t'.join(all_group_id) + '\n')
for each_1 in all_group_id:
    row = [each_1]
    for each_2 in all_group_id:
        current_key = '%s,%s' % (each_2, each_1)
        if current_key not in transfer_count:
            row.append('0')
        else:
            row.append(transfer_count[current_key])
    matrix_file.write('\t'.join(row) + '\n')
matrix_file.close()

os.remove('tmp1.txt')
os.remove('tmp1_sorted.txt')
os.remove('tmp1_sorted_count.txt')
print('Done!')
