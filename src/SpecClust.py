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
    sgaAsga2aff[(sga1, sga2)] = aff
    set_gene.add(sga1)
    set_gene.add(sga2)
    #print >> fo, l
    
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
        W[i][j] = sgaAsga2aff[(id2gene[i], id2gene[j])]


##
'''
W = np.array([[1,2,2,0,0,0],
    [2,1,2,0,0,0],
    [2,2,1,0,0,0],
    [0,0,0,1,8,0],
    [0,0,0,8,6,0],
    [0,0,0,0,0,1]])
'''
D = W.sum(axis=0)
Dinv = np.diag(1.0/D)

D = np.diag(D);
L =  D - W
#L = D - W;
L = np.dot(Dinv,L)

lambda0, V0 = linalg.eig(L)
#L = diag(1./sum(W))*L;
#[V0, lambda0] = eig(L);
ordr =lambda0.argsort()
lam = lambda0[ordr]
#V = V0
V = V0[:,ordr]
plt.plot(lam)
plt.ylabel('labmda')
plt.xlabel('k')
plt.title('SGA in BRCA')
plt.ylim((0,1.5))
#lam,ordr = sort(lambda0);
#lambda = diag(lam);
#v = V0(:,ord);


print 'Done!'
#fo.close()