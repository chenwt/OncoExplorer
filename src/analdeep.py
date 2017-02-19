# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 21:49:41 2017

@author: Yifeng Tao
"""

# Visualize and plot based on deepwalk.

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

# Generate enrichment file.

                        #nw = 10,20,50,100,200
# wl = 1,2,3,5,10,20,40

maxIter = 10
NW = [10,20,50,100,200]
#NW = [200]
WL = [1,2,3,5,10,20,40]
#WL = [1]
len_nw, len_wl = len(NW), len(WL)

enrich_true = np.zeros([len_wl,len_nw], dtype=float)
enrich_fake = np.zeros([len_wl,len_nw], dtype=float)

gene2goid = dd(set)
path = 'gene2goId.txt'
f = open(path, 'r')
#next(f)
for line in f:
    l = line.strip().split('\t')
    g, goid = l[0], l[1]
    gene2goid[g].add(goid)
f.close()

gid2gene = dd()
f = open('LUT_DW_brca.txt', 'r')
rg_gid = 0
for line in f:
    
    l = line.strip().split('\t')
    gene, gid = l[0], int(l[1])
    
    if gene.endswith('_deg'): break
    rg_gid += 1
    gid2gene[rg_gid] = gene
    #print gene, gid, rg_gid
f.close()

W = np.zeros([rg_gid, 40], float)


for count_wl,wl in enumerate(WL):
    #print count_wl
    for count_nw,nw in enumerate(NW):
        #TODO: note the number of id...
        path = '../result/brca_km40_wl'+str(wl)+'_nw'+str(nw)+'.embeddings'
        f = open(path, 'r')
        next(f)
        k = 0
        for line in f:
            l = line.strip().split(' ', 1)
            gid = int(l[0])
            if gid <= rg_gid:
                k += 1
                x = [float(i) for i in l[1].split(' ')]
                W[gid-1,:] = x
                #print k, gid, x
        f.close()
        
        X = W
        km = 20
        
        for _ in range(maxIter):       
            
                        
            
            kmeans = KMeans(n_clusters=km).fit(X)
            y = kmeans.labels_

                        
            
            
            
            y_fake = copy.deepcopy(y)
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
                    if y[i] == y[j]:
                        within_d += 1
                        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                            within_u += 1
            #print within_d, total_d, within_u, total_u
            print 'true', wl, nw, 1.0*within_u/within_d*total_d/total_u
            enrich_true[count_wl,count_nw] += 1.0*within_u/within_d*total_d/total_u
                    
            
            
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
            print 'fake', wl, nw, 1.0*within_u/within_d*total_d/total_u
            enrich_fake[count_wl,count_nw] += 1.0*within_u/within_d*total_d/total_u
    #enrich_ctrl[t] += 1.0*within_u/within_d*total_d/total_u
    #t += 1


                        #nw = 10,20,50,100,200
# wl = 1,2,3,5,10,20,40
np.save('true_km20'+'_brca',1.0*enrich_true/maxIter)
np.save('fake_km20'+'_brca',1.0*enrich_fake/maxIter)





'''


iter = 1
NW = [10,20,50,100,200]
WL = [1,2,3,5,10,20,40]
len_nw, len_wl = len(NW), len(WL)

enrich_true = np.zeros([WL,len_nw], dtype=float)
enrich_fake = np.zeros([WL,len_nw], dtype=float)

gene2goid = dd(set)
path = 'gene2goId.txt'
f = open(path, 'r')
#next(f)
for line in f:
    l = line.strip().split('\t')
    g, goid = l[0], l[1]
    gene2goid[g].add(goid)
f.close()

gid2gene = dd()
f = open('LUT_DW_brca.txt', 'r')
rg_gid = 0
for line in f:
    
    l = line.strip().split('\t')
    gene, gid = l[0], int(l[1])
    
    if gene.endswith('_deg'): break
    rg_gid += 1
    gid2gene[rg_gid] = gene
    #print gene, gid, rg_gid
f.close()

W = np.zeros([rg_gid, 40], float)

count_wl = 0
for wl in WL:
    count_wl += 1
    count_nw = 0
    for nw in NW:
        count_nw += 1

    #TODO: note the number of id...
    path = '../result/brca_km40_wl40_nw200.embeddings'
    f = open(path, 'r')
    next(f)
    k = 0
    for line in f:
        l = line.strip().split(' ', 1)
        gid = int(l[0])
        if gid <= rg_gid:
            k += 1
            x = [float(i) for i in l[1].split(' ')]
            W[gid-1,:] = x
            #print k, gid, x
    f.close()
    
    X = W
    km = 40
    kmeans = KMeans(n_clusters=km).fit(X)
    y = kmeans.labels_
    
    y_fake = copy.deepcopy(y)
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
            if y[i] == y[j]:
                within_d += 1
                if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                    within_u += 1
    #print within_d, total_d, within_u, total_u
    print 'true', km, 1.0*within_u/within_d*total_d/total_u
    
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
    print 'fake', km, 1.0*within_u/within_d*total_d/total_u



model = TSNE(n_components=2, random_state=0)
#np.set_printoptions(suppress=True)
Xnew = model.fit_transform(X)
h = Xnew[:,0]
v = Xnew[:,1]

#for i in range(rg_gid):
#    sga1 = set_gene[i]
#    for line in sga2sgaAaff[sga1]:
#        sga2, aff = line[0], line[1]
##        sga1id = gene2id[sga1]
##        sga2id = gene2id[sga2]
##        gene2degree[sga1] += aff
##        gene2degree[sga2] += aff
#        ax.plot([h[i], h[sga2id]], [v[sga1id], v[sga2id]], marker='.', linestyle='-', color = [0.7, 0.7, 0.7] )

RGB_tuples = _get_colors(km)
fig, ax = plt.subplots(figsize=(8, 8))

for i in range(rg_gid):
    ax.plot(h[i], v[i], marker='o', linestyle='', markersize = 5, color = RGB_tuples[y[i]])

'''
print 'Done!'





