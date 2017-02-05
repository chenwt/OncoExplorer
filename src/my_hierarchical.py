# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd
import numpy as np
#import scipy.io
import os
from scipy.io import savemat
from scipy.sparse import csr_matrix
import time

from scipy.cluster.hierarchy import linkage
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

#def augmented_dendrogram(*args, **kwargs):
#
#    ddata = dendrogram(*args, **kwargs)
#
#    #if not kwargs.get('no_plot', False):
#    for i, d in zip(ddata['icoord'], ddata['dcoord']):
#        x = 0.5 * sum(i[1:3])
#        y = d[1]
#        plt.plot(x, y, 'ro')
#        plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
#                     textcoords='offset points',
#                     va='top', ha='center')
#
#    return ddata
    
# Generate a random sample of `n` points in 2-d.
    
    
    
sga_label = set()
sgaAsga2dist = dd(float)
xx2dist = dd(float)
path = 'brca_baseline1.txt'
f = open(path, 'r')
for line in f:
    line = line.strip().split('\t')
    sga1, sga2, aff = line[0], line[1], float(line[2])
    sga_label.add(sga1)
    sga_label.add(sga2)
    sgaAsga2dist[(sga1,sga2)] = -aff
    #print sga1, sga2, aff
f.close()



sga_label = list(sga_label)
sga2x = dd(float)
x2sga = dd(set)
x = 0.0
for sga in sga_label:
    sga2x[sga] = x
    x2sga[x].add(sga)
    x += 1

for line in sgaAsga2dist.keys():
    sga1,sga2,dist = line[0], line[1], sgaAsga2dist[line]
    xx2dist[(sga2x[sga1], sga2x[sga2])] = dist

# Merging algorithm.
# S1: Set that maintain the remaining sets
set_x = set()
for i in range(len(sga_label)):
    set_x.add(i)

#f = open('out.txt', 'w')

#def distBetweenSets(S1, S2):
#    # max mode
#    dist = -1000.0
#    for s1 in S1:
#        for s2 in S2:
#            if s1 == s2: print 'error!', s1, s2
#            tmp = sgaAsga2dist[(s1, s2)]
#            if tmp > dist: dist = tmp
#    return dist

#def distBetweenSets(S1, S2):
#    # min mode
#    dist = 0
#    for s1 in S1:
#        for s2 in S2:
#            if s1 == s2: print 'error!', s1, s2
#            tmp = sgaAsga2dist[(s1, s2)]
#            if tmp < dist: dist = tmp
#    return dist

def distBetweenSets(S1, S2):
    # median mode
    
    dist = 0
    #l1, l2 = len(S1), len(S2)
    m = 10
    dist_list = np.ones(m, dtype=float)
    
    for s1 in S1:
        for s2 in S2:
            if s1 == s2: print 'error!', s1, s2
            tmp = sgaAsga2dist[(s1, s2)]
            for k in range(m):
                if tmp < dist_list[k]:
                    for t in range(m-1, k, -1):
                        dist_list[t] = dist_list[t-1]
                    dist_list[k] = tmp
                    break
    for t in range(m-1, -1, -1):
        dist = dist_list[t]
        if dist <= 0: break
        
    #dist = dist_list[m-1]
    #print dist_list, dist
            #dist_list[k] = tmp
    #dist_list = np.array(dist_list)
    #dist = np.median(dist_list)
    
        #if tmp < dist: dist = tmp
    return dist

time_start = time.time()
link_matrix = np.zeros((len(sga_label)-1, 4), dtype=float)
#for count in range(1):
for count in range(len(sga_label)-1):
    # introduce random
    min_dist = 1.0
    min_x1 = -1.0
    min_x2 = -1.0
    for x1 in set_x:
        #print x1
        for x2 in set_x:
            if x1 == x2: continue
            S1, S2 = x2sga[x1], x2sga[x2]
            tmp = distBetweenSets(S1, S2)
            if tmp < min_dist:
                min_x1 = x1
                min_x2 = x2
                min_dist = tmp
    set_x.remove(min_x1)
    set_x.remove(min_x2)
    xnew = len(sga_label)+count
    set_x.add(xnew)
    newS = x2sga[min_x1].union(x2sga[min_x2])
    x2sga[xnew] = newS
    link_matrix[count][0] = min_x1
    link_matrix[count][1] = min_x2
    link_matrix[count][2] = min_dist
    link_matrix[count][3] = len(newS)
    print min_x1, min_x2, min_dist, xnew
    print x2sga[min_x1], x2sga[min_x2]
            #print x1, x2, dist
            #print >> f, str(x1)
    #pass
#f.close()




time_end = time.time()
print time_end-time_start
for count in range(len(sga_label)-1):
    link_matrix[count][2] += 205
#print 
#linkage_matrix
plt.figure(2, figsize=(4, 50))
plt.clf()

plt.subplot(1, 2, 1)
ddata = dendrogram(link_matrix, color_threshold=1, 
                   show_leaf_counts=True,
                   labels=sga_label, orientation = 'left', leaf_font_size = 10)

'''
    
np.random.seed(12312)
n = 100
x = np.random.multivariate_normal([0, 0], np.array([[4.0, 2], [2, 1.4]]),
                                  size=(n,))

plt.figure(1, figsize=(6, 5))
plt.clf()
plt.scatter(x[:, 0], x[:, 1])
plt.axis('equal')
plt.grid(True)


linkage_matrix = linkage(x, "single")

plt.figure(2, figsize=(4, 50))
plt.clf()

plt.subplot(1, 2, 1)
show_leaf_counts = False

lb = []
for i in range(n):
    lb.append('a'+str(i))

ddata = dendrogram(linkage_matrix, color_threshold=1, 
                   show_leaf_counts=show_leaf_counts,
                   labels=lb, orientation = 'left', leaf_font_size = 10)
#ddata = dendrogram(linkage_matrix, color_threshold=1, p=6,
#               truncate_mode='lastp', show_leaf_counts=show_leaf_counts,
#               labels=lb)
#ddata = augmented_dendrogram(linkage_matrix,
#               color_threshold=1,
#               p=6,
#               truncate_mode='lastp',
#               show_leaf_counts=show_leaf_counts,
#               )
               
plt.title("show_leaf_counts = %s" % show_leaf_counts)

#plt.subplot(1, 2, 2)
#show_leaf_counts = True
#ddata = augmented_dendrogram(linkage_matrix,
#               color_threshold=1,
#               p=6,
#               truncate_mode='lastp',
#               show_leaf_counts=show_leaf_counts,
#               )
#plt.title("show_leaf_counts = %s" % show_leaf_counts)

plt.show()

'''
#EOF.
