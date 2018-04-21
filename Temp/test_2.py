import os
from datetime import datetime
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from ete3 import Tree, ClusterTree



wd = '/Users/songweizhi/Desktop/tree_135_138'

file_in = '/Users/songweizhi/Desktop/tree_135_138/taxon_assignment_lineage_c80.txt'



rank = 'd'


def get_rank_assignment(rank, taxon_assignment_lineage_file):

    rank_to_position_dict = {'d': 1, 'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6, 's': 7}
    rank_position = rank_to_position_dict[rank]

    assignment_dict = {}
    for each in open(taxon_assignment_lineage_file):
        each_split = each.strip().split('\t')
        bin_name = each_split[0]
        assignment = each_split[1].split(';')
        assignment_no_num = []
        for each_assign in assignment:
            assignment_no_num.append(each_assign.split('(')[0])

        new_assign = ''
        if len(assignment_no_num) <= rank_position:
            new_assign = assignment_no_num[-1]
        elif len(assignment_no_num) > rank_position:
            new_assign = assignment_no_num[rank_position-1]

        assignment_dict[bin_name] = new_assign

    return assignment_dict





print(get_rank_assignment(rank, file_in))