# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 04:52:07 2017

@author: Yifeng Tao
"""

# Generate kNN graph.



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

itertime = 20
KM = np.array([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
#KM = np.array([40])
classsize = np.zeros([0,], dtype=float)
enrich = np.zeros([len(KM),],dtype = float)
enrich_ctrl = np.zeros([len(KM),],dtype = float)

def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


#for iter in range(itertime):
kNN = 5
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
    #sgaAsga2aff[(sga1,sga2)] = aff        
    
    set_gene.add(sga1)
    set_gene.add(sga2)
f.close()

set_gene = list(set_gene)

for sga in set_gene:
    context = sga2sgaAaff[sga]
    context = list(context)    
    # shuffle
    
    for i in range(kNN):
        shuffle(context)
    
        Maff = 0.0
        Mgene = 'helloword'
        for line in context:
            sgacontext, aff = line[0], line[1]
            if aff > Maff:
                Maff = aff
                Mgene = sgacontext
        sgaAsga2aff[(sga, Mgene)] = aff
        sgaAsga2aff[(Mgene, sga)] = aff
        #print sga, Mgene, Maff
        context.remove((Mgene, Maff))
        if len(context) == 0: break    



#set_gene = list(set_gene)
    
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


path = 'sga2sgaAff_baseline_kNN_'+str(kNN)+'_brca.txt'
f = open(path, 'w')
print >> f, 'Source\tTarget\tWeight'
for line in sgaAsga2aff.keys():
        sga1, sga2 = line[0], line[1]
        overlap = math.sqrt(1.0*sgaAsga2aff[(sga1, sga2)]*sgaAsga2aff[(sga2, sga1)])
        #sga2, overlap = line[0], line[1]
        if overlap == 0: continue
        if sga1 == sga2: continue
        print >> f, sga1 + '\t' + sga2 + '\t' + str(overlap)
        #print >> f, sga1 + '\t' + sga2 + '\t' + overlap
f.close()


sga2sgaAaff = dd(set)
set_gene = set()

path = 'sga2sgaAff_baseline_kNN_'+str(kNN)+'_brca.txt'
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
#L = D - W;
L = np.dot(Dinv,L)

#print 2

lambda0, V0 = linalg.eig(L)

#print 3
#L = diag(1./sum(W))*L;
#[V0, lambda0] = eig(L);
ordr =lambda0.argsort()
lam = lambda0[ordr]
#V = V0
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
    
    
  






    
km = 40
#km = 50
X = V[:,0:km]

sga2sgaAaff = dd(set)

W = X
M = 0
for id1 in range(len_gene):
    #g1 = id2gene[id1]
    for id2 in range(len_gene):
        #g2 = id2gene[id2]
        d = W[id1,:] - W[id2,:]
        d =  np.linalg.norm(d)
        M = max(d, M)
        
for id1 in range(len_gene):
    g1 = id2gene[id1]
    for id2 in range(len_gene):
        g2 = id2gene[id2]
        d = W[id1,:] - W[id2,:]
        d =  M - np.linalg.norm(d)
        
        sga2sgaAaff[g1].add((g2, d))





def getsgaAsga2pp():
    set_pathway = set()
    sgaAsga2pp = dd(float)
    path = '../src/pp/sgaAsga2pp.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sga1, sga2, prob = l[0], l[1], float(l[2])
        sgaAsga2pp[(sga1, sga2)] = prob
        set_pathway.add(sga1)
        set_pathway.add(sga2)
    
    f.close()
    return sgaAsga2pp, set_pathway

sgaAsga2pp, set_pathway = getsgaAsga2pp()

avgacc = 0
for trial in range(10):
    count_total = 0
    count_true = 0
    for sga in set_gene:
        if sga not in sga2sgaAaff.keys(): continue
        count_total += 1
        NNgene = 'helloworld'
        NNdist = 0.0
        Context = list(sga2sgaAaff[sga])
        shuffle(Context)
        for line in Context:
            sgaContext, aff = line[0], line[1]
            if sgaContext == sga: continue
            if aff > NNdist:
                NNdist = aff
                NNgene = sgaContext
        #print sga, NNgene, NNdist
        thd = 0.001
        if sgaAsga2pp[(sga, NNgene)] == 0.0: continue
        if sgaAsga2pp[(sga, NNgene)] <= thd:
            count_true += 1
#        ovlp = gene2goid[sga].intersection(gene2goid[NNgene])
#        ovlp = len(ovlp)
#        metric = 'one'
#        if metric == 'one':
#            if ovlp > 0:
#                count_true += 1
#        elif metric == 'total':
#            count_true += ovlp
#        elif metric == 'jaccard':
#            count_true += 1.0*ovlp/len(gene2goid[sga].union(gene2goid[NNgene]))
    acc = 1.0*count_true/count_total
    avgacc += acc
    # Omit if necessary
    #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
print 'cancer:brca\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)


print 'Done!'


#EOF.