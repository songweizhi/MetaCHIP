
# needed imports
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
import numpy as np



def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

# some setting for this notebook to actually show the graphs inline
# you probably won't need this
np.set_printoptions(precision=5, suppress=True)  # suppress scientific float notation


# generate two clusters: a with 100 points, b with 50:
np.random.seed(4711)  # for repeatability of this tutorial
a = np.random.multivariate_normal([10, 0], [[3, 1], [1, 4]], size=[100,])
b = np.random.multivariate_normal([0, 20], [[3, 1], [1, 4]], size=[50,])
X = np.concatenate((a, b),)

print(X.shape)  # 150 samples with 2 dimensions
print(type(X))
plt.scatter(X[:,0], X[:,1])
#plt.show()


# generate the linkage matrix
Z = linkage(X, 'ward')
c, coph_dists = cophenet(Z, pdist(X))

print('c')
print(c)
print('coph_dists')
print(coph_dists)
print(type(coph_dists))


# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(Z, leaf_rotation=90, leaf_font_size=8)
plt.show()


plt.title('Hierarchical Clustering Dendrogram (truncated)')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(Z,truncate_mode='lastp',  # show only the last p merged clusters
    p=12,  # show only the last p merged clusters
    show_leaf_counts=False,  # otherwise numbers in brackets are counts
    leaf_rotation=90.,leaf_font_size=12.,show_contracted=True,  # to get a distribution impression in truncated branches
)
plt.show()


plt.title('Hierarchical Clustering Dendrogram (truncated)')
plt.xlabel('sample index or (cluster size)')
plt.ylabel('distance')
dendrogram(Z,truncate_mode='lastp',  # show only the last p merged clusters
    p=12,  # show only the last p merged clusters
    leaf_rotation=90.,leaf_font_size=12.,show_contracted=True,  # to get a distribution impression in truncated branches
)
plt.show()


fancy_dendrogram(Z,truncate_mode='lastp',p=12,leaf_rotation=90.,leaf_font_size=12.,show_contracted=True,annotate_above=10,  # useful in small plots so annotations don't overlap
)
plt.show()

# set cut-off to 50
max_d = 50  # max_d as in max_distance

fancy_dendrogram(Z,truncate_mode='lastp',p=12,leaf_rotation=90.,leaf_font_size=12.,show_contracted=True,annotate_above=10,max_d=max_d,)
plt.show()

fancy_dendrogram(Z,truncate_mode='lastp',p=12,leaf_rotation=90.,leaf_font_size=12.,show_contracted=True, annotate_above=10,max_d=16,)
plt.show()


max_d = 23
clusters = fcluster(Z, max_d, criterion='distance')
print(clusters)


plt.figure(figsize=(10, 8))
plt.scatter(X[:,0], X[:,1], c=clusters, cmap='prism')  # plot points with cluster dependent colors
plt.show()


