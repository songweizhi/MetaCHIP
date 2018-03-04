import os
import argparse
import itertools
from string import ascii_uppercase

usage = """

python /srv/scratch/z5039045/Scripts/cluster_2_grouping_file.py -c grouping_out.txt -g grouping.txt

"""

def get_group_index_list():
    def iter_all_strings():
        size = 1
        while True:
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)
            size += 1

    group_index_list = []
    for s in iter_all_strings():
        group_index_list.append(s)
        if s == 'ZZ':
            break
    return group_index_list


def cluster_2_grouping_file(cluster_file, grouping_file):

    t1 = 'cluster_tmp1.txt'
    t1_sorted = 'cluster_tmp1_sorted.txt'

    t1_handle = open(t1, 'w')
    for each in open(cluster_file):
        if not each.startswith(','):
            genome_id = each.strip().split(',')[0]
            cluster_id = each.strip().split(',')[1]
            t1_handle.write('%s,%s\n' % (cluster_id, genome_id))
    t1_handle.close()

    os.system('cat %s | sort > %s' % (t1, t1_sorted))
    group_index_list = get_group_index_list()
    grouping_file_handle = open(grouping_file, 'w')

    current_cluster_name = ''
    group_index_no = 0
    n = 1
    for each in open(t1_sorted):
        cluster_name = each.strip().split(',')[0]
        genome_name = each.strip().split(',')[1]

        if current_cluster_name == '':
            current_cluster_name = cluster_name
            grouping_file_handle.write('%s_%s,%s\n' % (group_index_list[group_index_no], n, genome_name))
            n += 1
        elif current_cluster_name == cluster_name:
            grouping_file_handle.write('%s_%s,%s\n' % (group_index_list[group_index_no], n, genome_name))
            n += 1
        elif current_cluster_name != cluster_name:
            current_cluster_name = cluster_name
            group_index_no += 1
            n = 1
            grouping_file_handle.write('%s_%s,%s\n' % (group_index_list[group_index_no], n, genome_name))
            n += 1

    os.remove(t1)
    os.remove(t1_sorted)


parser = argparse.ArgumentParser()

parser.add_argument('-c',
                    required=True,
                    help='cluster file')

parser.add_argument('-g',
                    required=True,
                    help='grouping file')

args = vars(parser.parse_args())

cluster_file = args['c']
grouping_file = args['g']

cluster_2_grouping_file(cluster_file, grouping_file)

