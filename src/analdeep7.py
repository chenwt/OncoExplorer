# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 04:39:27 2017

@author: Yifeng Tao
"""

#analdeep3_1.py



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





gid2gene = dd()
gene2gid = dd()
f = open('LUT_DW_brca.txt', 'r')
rg_gid = 0
for line in f:
    
    l = line.strip().split('\t')
    gene, gid = l[0], int(l[1])
    
    if gene.endswith('_deg'): break
    rg_gid += 1
    gid2gene[rg_gid-1] = gene
    gene2gid[gene] = rg_gid-1
    #print gene, gid, rg_gid
f.close()

Freq = [2,3,4,6,8,10, 20, 30, 40, 50]

for i in range(30):
    Freq.append(10*i + 60)

l_g = len(gid2gene)
lf = len(Freq)
treat = np.zeros([lf,])
ctrl = np.zeros([lf,])
Iter = 5
#for freq in [3, 5, 10, 20, 40, 80, 200]:
for it, freq in enumerate(Freq):

    y = np.zeros([l_g,])
    y = [(i % freq) for i in range(l_g)]
    gene2goid = dd(set)
    for i in range(l_g):
        gene2goid[gid2gene[i]].add(i%freq)
        
    
    for itera in range(Iter):
        
        total_u = 0.0
        total_d = 0.0
        within_u = 0.0
        within_d = 0.0
        for i in range(rg_gid):
            g1 = gid2gene[i]
            for j in range(rg_gid):
                if i == j: continue
                g2 = gid2gene[j]
                total_d += 1
                if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                    total_u += 1
                if y[i] == y[j]:
                    within_d += 1
                    if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                        within_u += 1
        
        enrich =1.0*within_u/within_d*total_d/total_u
        if itera == 0:
            treat[it] = enrich
        else:
            ctrl[it] += enrich               
        #print within_d, total_d, within_u, total_u
        print '#cluster=',freq, 'enrichment=', enrich
        #print 'iter=',itera,'size=',529/freq, 'enrichment=', 1.0*within_u/within_d*total_d/total_u
        #enrich_true[count_wl,count_nw] += 1.0*within_u/within_d*total_d/total_u   
        shuffle(y)
    #W = np.zeros([rg_gid, 40], float)
ctrl = 1.0*ctrl/(Iter-1)

#plt.plot(Freq, treat, linestyle='-', linewidth=2,color=RGB_tuples[count_wl], label=WL[count_wl])
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

csize = [i for i in Freq]
ax.plot(csize, treat, linestyle='-', linewidth=2,marker='o', markersize=10,)
ax.plot(csize, ctrl, linestyle='-', linewidth=2,marker='o', markersize=10,)
ax.set_yscale('log')
plt.xlabel('cluster size')
plt.ylabel('enrichment')
plt.legend(['true cluster', 'randomly shuffled'])
plt.figure()
plt.plot(csize, 1.0*treat/ctrl, linestyle='-', linewidth=2,marker='o', markersize=10,)
plt.xlabel('cluster size')
plt.ylabel('enrichment ratio')

#plt.title('Enrichment of GoTerm (Deep walk), 40-means x10')
#plt.xlabel('#walks/#nodes')
#plt.ylabel('enrichment')
#plt.xlim([0.0,210])
#plt.ylim([1.0,2.0])
#plt.legend(le, WL,loc=3)

print 'Done!'


#EOF.