import io
import os
import argparse
import itertools
import pandas as pd
from ete3 import Tree
from string import ascii_uppercase
from scipy.cluster.hierarchy import linkage
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np


def get_node_distance(tree, node_1, node_2):
    distance = tree.get_distance(node_1, node_2)
    return distance


def get_matrix(input_newick_file, output_matrix_file):
    # read in tree
    tree_in = Tree(input_newick_file, format=3)

    # get sorted leaf node list
    leaf_node_list = []
    for leaf_node in tree_in:
        leaf_node_list.append(leaf_node.name)

    leaf_node_list = sorted(leaf_node_list)

    # write out distance matrix
    matrix_out_handle = open(output_matrix_file, 'w')
    matrix_out_handle.write(',%s\n' % ','.join(leaf_node_list))
    for each_node_1 in leaf_node_list:
        current_node_distance_list = []
        for each_node_2 in leaf_node_list:
            distance = 0
            if each_node_1 != each_node_2:
                distance = get_node_distance(tree_in, each_node_1, each_node_2)
                distance = float("{0:.5f}".format(distance))
            current_node_distance_list.append(str(distance))
        matrix_out_handle.write('%s,%s\n' % (each_node_1, ','.join(current_node_distance_list)))

    matrix_out_handle.close()

os.chdir('/Users/songweizhi/Desktop')

newick_tree_file = 'species_tree.newick'
matrix_file = 'species_tree_with_group_matrix.csv'

# get distance matrix from tree
get_matrix(newick_tree_file, matrix_file)

# read in distance matrix as data frame
distance_matrix = pd.read_csv(matrix_file, index_col=0)

# read in distance matrix as numpy array
# import numpy as np
# distance_array = np.genfromtxt(matrix_file, delimiter=',', skip_header=True)
#
# print(distance_array.shape)
# # run clustering
#

#
# from scipy.cluster.hierarchy import dendrogram, linkage
# from matplotlib import pyplot as plt
#
#
# np.random.seed(4711)  # for repeatability of this tutorial
# a = np.random.multivariate_normal([10, 0], [[3, 1], [1, 4]], size=[100,])
# b = np.random.multivariate_normal([0, 20], [[3, 1], [1, 4]], size=[50,])
# X = np.concatenate((a, b),)
#
#
# print(a)
#
# print(X.shape)
# print(type(X))
#
# print(distance_matrix.shape)
# print(type(distance_matrix))



# Z = linkage(distance_matrix)
# Z
# dendrogram(Z)
#
# plt.show()


#
# def kMedoids(D, k, tmax=100):
#     # determine dimensions of distance matrix D
#     m, n = D.shape
#
#     if k > n:
#         raise Exception('too many medoids')
#
#     # find a set of valid initial cluster medoid indices since we
#     # can't seed different clusters with two points at the same location
#     valid_medoid_inds = set(range(n))
#     invalid_medoid_inds = set([])
#     rs,cs = np.where(D==0)
#     # the rows, cols must be shuffled because we will keep the first duplicate below
#     index_shuf = range(len(rs))
#     np.random.shuffle(index_shuf)
#     rs = rs[index_shuf]
#     cs = cs[index_shuf]
#     for r,c in zip(rs,cs):
#         # if there are two points with a distance of 0...
#         # keep the first one for cluster init
#         if r < c and r not in invalid_medoid_inds:
#             invalid_medoid_inds.add(c)
#     valid_medoid_inds = list(valid_medoid_inds - invalid_medoid_inds)
#
#     if k > len(valid_medoid_inds):
#         raise Exception('too many medoids (after removing {} duplicate points)'.format(
#             len(invalid_medoid_inds)))
#
#     # randomly initialize an array of k medoid indices
#     M = np.array(valid_medoid_inds)
#     np.random.shuffle(M)
#     M = np.sort(M[:k])
#
#     # create a copy of the array of medoid indices
#     Mnew = np.copy(M)
#
#     # initialize a dictionary to represent clusters
#     C = {}
#     for t in xrange(tmax):
#         # determine clusters, i. e. arrays of data indices
#         J = np.argmin(D[:,M], axis=1)
#         for kappa in range(k):
#             C[kappa] = np.where(J==kappa)[0]
#         # update cluster medoids
#         for kappa in range(k):
#             J = np.mean(D[np.ix_(C[kappa],C[kappa])],axis=1)
#             j = np.argmin(J)
#             Mnew[kappa] = C[kappa][j]
#         np.sort(Mnew)
#         # check for convergence
#         if np.array_equal(M, Mnew):
#             break
#         M = np.copy(Mnew)
#     else:
#         # final update of cluster memberships
#         J = np.argmin(D[:,M], axis=1)
#         for kappa in range(k):
#             C[kappa] = np.where(J==kappa)[0]
#
#     # return results
#     return M, C
#
#
#
# M,C = kMedoids(distance_matrix, 10, tmax=100)
#
# print(M)
# print(C)
