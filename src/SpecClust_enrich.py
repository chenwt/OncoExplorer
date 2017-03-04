# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:48:20 2017

@author: Yifeng Tao
"""

# Spectral clustering.

import numpy as np
import numpy.linalg as linalg
from collections import defaultdict as dd
import matplotlib.pyplot as plt
from random import shuffle
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.manifold import TSNE
import colorsys
import math
import copy



def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


    
sga2sgaAaff = dd(set)
sgaAsga2aff = dd(int)
set_gene = set()
path = 'sga2sgaAff_baseline_brca.txt'
f = open(path, 'r')
pathout = 'out.txt'
fo = open(pathout, 'w')
next(f)
for line in f:
    l = line.strip().split('\t')
    sga1, sga2, aff = l[0], l[1], int(l[2])
    sga2sgaAaff[sga1].add((sga2, aff))
    sgaAsga2aff[(sga1,sga2)] = aff        
    
    set_gene.add(sga1)
    set_gene.add(sga2)
f.close()

set_gene = list(set_gene)
    
gene2id = dd()
id2gene = dd()

k = 0
for g in set_gene:
    gene2id[g] = k
    id2gene[k] = g
    k += 1


len_gene = len(set_gene)
W = np.zeros([len_gene, len_gene], float)
for i in range(len_gene):
    for j in range(len_gene):
        W[i][j] = math.sqrt(1.0*sgaAsga2aff[(id2gene[i], id2gene[j])]*sgaAsga2aff[(id2gene[j], id2gene[i])])


path = 'sga2sgaAff_baseline_simple_brca.txt'
f = open(path, 'w')
print >> f, 'Source\tTarget\tWeight'
for line in sgaAsga2aff.keys():
        sga1, sga2 = line[0], line[1]
        overlap = math.sqrt(1.0*sgaAsga2aff[(sga1, sga2)]*sgaAsga2aff[(sga2, sga1)])
        if overlap == 0: continue
        if sga1 == sga2: continue
        print >> f, sga1 + '\t' + sga2 + '\t' + str(overlap)
f.close()


sga2sgaAaff = dd(set)
set_gene = set()

path = 'sga2sgaAff_baseline_simple_brca.txt'
f = open(path, 'r')
next(f)
for line in f:
    l = line.strip().split('\t')
    sga1, sga2, aff = l[0], l[1], float(l[2])
    sga2sgaAaff[sga1].add((sga2, aff))
    set_gene.add(sga1)
    set_gene.add(sga2)
f.close()

set_gene = list(set_gene)

gene2id = dd()
id2gene = dd()

k = 0
for g in set_gene:
    gene2id[g] = k
    id2gene[k] = g
    k += 1


len_gene = len(set_gene)
W = np.zeros([len_gene, len_gene], float)
for i in range(len_gene):
    for j in range(len_gene):
        W[i][j] = math.sqrt(1.0*sgaAsga2aff[(id2gene[i], id2gene[j])]*sgaAsga2aff[(id2gene[j], id2gene[i])])

D = W.sum(axis=0)
Dinv = np.diag(1.0/D)

D = np.diag(D);
L =  D - W
L = np.dot(Dinv,L)

lambda0, V0 = linalg.eig(L)

ordr =lambda0.argsort()
lam = lambda0[ordr]
V = V0[:,ordr]



gene2goid = dd(set)
path = 'gene2goId.txt'
f = open(path, 'r')
#next(f)
for line in f:
    l = line.strip().split('\t')
    g, goid = l[0], l[1]
    gene2goid[g].add(goid)
f.close()




t = 0
km = 40
X = V[:,0:km]

kmeans = KMeans(n_clusters=km).fit(X)
y = kmeans.labels_

c2gene = dd(set)
path = '../src/pp/cluster.txt'
f = open(path, 'w')
for i in range(len_gene):
    genet = id2gene[i]
    print >> f, genet+'\t'+str(y[i])
    c2gene[y[i]].add(genet)
f.close()

path = '../src/pp/c2gene.txt'
f = open(path, 'w')
for i in range(km):
    #print >> f, i
    for genet in c2gene[i]:
        print >> f, str(i)+'\t'+genet
f.close()

path = '../src/pp/c2gene1.txt'
f = open(path, 'w')
for i in range(km):
    print >> f, i
    for genet in c2gene[i]:
        print >> f, genet
f.close()
#for i in range(len_gene):
#    for j in range(len_gene):
#        W[i][j] = math.sqrt(1.0*sgaAsga2aff[(, id2gene[j])]*sgaAsga2aff[(id2gene[j], id2gene[i])])
#







y_fake = copy.deepcopy(y)
shuffle(y_fake)

len_gene = len(set_gene)

total_u = 0.0
total_d = 0.0
within_u = 0.0
within_d = 0.0
for i in range(len_gene):
    g1 = id2gene[i]
    for j in range(len_gene):
        if i == j: continue
        g2 = id2gene[j]
        l1 = y[i]
        l2 = y[j]
        total_d += 1
        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
            total_u += 1
        if y[i] == y[j]:
            within_d += 1
            if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                within_u += 1
print 'true', km, 1.0*within_u/within_d*total_d/total_u


total_u = 0.0
total_d = 0.0
within_u = 0.0
within_d = 0.0
for i in range(len_gene):
    g1 = id2gene[i]
    for j in range(len_gene):
        if i == j: continue
        g2 = id2gene[j]

        total_d += 1
        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
            total_u += 1
        if y_fake[i] == y_fake[j]:
            within_d += 1
            if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                within_u += 1
#print within_d, total_d, within_u, total_u
print 'fake', km, 1.0*within_u/within_d*total_d/total_u



model = TSNE(n_components=2)
#np.set_printoptions(suppress=True)
Xnew = model.fit_transform(X)
h = Xnew[:,0]
v = Xnew[:,1]




fig, ax = plt.subplots(figsize=(8, 8))




#set_gene = list(set_gene)



gene2degree = dd(float)
for i in range(len(set_gene)):
    sga1 = set_gene[i]
    for line in sga2sgaAaff[sga1]:
        sga2, aff = line[0], line[1]
        sga1id = gene2id[sga1]
        sga2id = gene2id[sga2]
        gene2degree[sga1] += aff
        gene2degree[sga2] += aff
        ax.plot([h[sga1id], h[sga2id]], [v[sga1id], v[sga2id]], marker='.', linestyle='-', color = [0.7, 0.7, 0.7] )

Mval = 0
for g in gene2degree.keys():
    degree = gene2degree[g]
    if degree > Mval: Mval = degree

Mval = math.sqrt(Mval)
for g in gene2degree.keys():
    degree = gene2degree[g]
    gene2degree[g] = math.sqrt(degree)/Mval



RGB_tuples = _get_colors(km)


for i in range(len(set_gene)):
    sga1 = set_gene[i]
    sga1id = gene2id[sga1]
    ax.plot(h[sga1id], v[sga1id], marker='o', linestyle='', markersize = int(20*gene2degree[sga1])+5, color = RGB_tuples[y[sga1id]])

plt.title('t-SNE of spectral clustering (k = '+str(km)+')')


print 'Done!'

#EOF.