# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 09:34:03 2017

@author: Yifeng Tao
"""

# Simulate and test the enrichment of random shuffle.

# We check result of DW, since it has high random enrichment.

import numpy as np
from collections import defaultdict as dd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

gid2gene = dd()
f = open('LUT_DW_brca.txt', 'r')
len_g = 0
for line in f:
    l = line.strip().split('\t')
    gene, gid = l[0], int(l[1])
    if gene.endswith('_deg'): break
    len_g += 1
    gid2gene[len_g] = gene
f.close()

#print len_g #529

gene2goid = dd(set)
path = 'gene2goId.txt'
f = open(path, 'r')
#next(f)
for line in f:
    l = line.strip().split('\t')
    g, goid = l[0], l[1]
    gene2goid[g].add(goid)
f.close()

wl = 1
nw = 200
path = '../result2/brca_km40_ss_wl'+str(wl)+'_nw'+str(nw)+'.embeddings'
f = open(path, 'r')
for line in f:
    l = line.strip().split(' ')
    len_v = int(l[1])
    break

W = np.zeros([len_g, len_v], dtype='float')

for line in f:
    l = line.strip().split(' ', 1)
    gid = int(l[0])
    #rg_gid =  len(l[1].split(' '))
    if gid <= len_g:
        x = [float(i) for i in l[1].split(' ')]
        W[gid-1,:] = x
        #print k, gid, x
f.close()

km = 40
kmeans = KMeans(n_clusters=km, random_state=0).fit(W)
y = kmeans.labels_

total_u = 0.0
total_d = 0.0
within_u = 0.0
within_d = 0.0
for i in range(len_g):
    g1 = gid2gene[i+1]
    for j in range(len_g):
        if i == j: continue # ORIGIN OF BUGs.
        g2 = gid2gene[j+1]
        total_d += 1
        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
            total_u += 1
        if y[i] == y[j]:
            within_d += 1
            if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                within_u += 1
#print within_d, total_d, within_u, total_u
within = 1.0*within_u/within_d
total = 1.0*total_u/total_d
print 'wl:',wl,'nw:',nw,'iter:',0,'true', within_d, within, total_d, total, 1.0*within/total




#enrich_true[count_wl,count_nw] += 1.0*within_u/within_d*total_d/total_u

'''
X = W
km = 40

#for itera in range(maxIter):



        
total_u = 0.0
total_d = 0.0
within_u = 0.0
within_d = 0.0
for i in range(rg_gid):
    g1 = gid2gene[i+1]
    for j in range(rg_gid):
        g2 = gid2gene[j+1]
        total_d += 1
        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
            total_u += 1
        if y[i] == y[j]:
            within_d += 1
            if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                within_u += 1
#print within_d, total_d, within_u, total_u
print 'wl:',wl,'nw:',nw,'iter:',itera,'true', wl, nw, 1.0*within_u/within_d*total_d/total_u
enrich_true[count_wl,count_nw] += 1.0*within_u/within_d*total_d/total_u
        
y_fake = y
shuffle(y_fake)           

total_u = 0.0
total_d = 0.0
within_u = 0.0
within_d = 0.0
for i in range(rg_gid):
    g1 = gid2gene[i+1]
    for j in range(rg_gid):
        g2 = gid2gene[j+1]
        total_d += 1
        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
            total_u += 1
        if y_fake[i] == y_fake[j]:
            within_d += 1
            if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                within_u += 1
#print within_d, total_d, within_u, total_u
print 'wl:',wl,'nw:',nw,'iter:',itera,'fake', wl, nw, 1.0*within_u/within_d*total_d/total_u
enrich_fake[count_wl,count_nw] += 1.0*within_u/within_d*total_d/total_u
'''

#EOF.