# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 04:40:34 2017

@author: Yifeng Tao
"""

# Analyze the result of spectral clustering.

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


meth = 3
RGB_tuples = _get_colors(meth)
#ratio = np.zeros((3,20),dtype=float)

plt.figure(figsize=(7,5))
le = list()

num_clus = np.linspace(5, 100, 20)
treat = np.load('enrichment_baseline_mukNN_'+str(10)+'_brca.npy')
plt.plot(num_clus, treat, marker='^', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
treat = np.load('enrichment_baseline_mukNN_'+str(20)+'_brca.npy')
plt.plot(num_clus, treat, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
treat = np.load('enrichment_baseline_mukNN_'+str(40)+'_brca.npy')
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
    #plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
#    markersize = 5, color = RGB_tuples[0])

treat = np.load('enrichment_baseline_kNN_'+str(3)+'_brca.npy')
plt.plot(num_clus, treat, marker='^', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
treat = np.load('enrichment_baseline_kNN_'+str(5)+'_brca.npy')
plt.plot(num_clus, treat, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
treat = np.load('enrichment_baseline_kNN_'+str(10)+'_brca.npy')
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
#plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
#    markersize = 5, color = RGB_tuples[1])

treat = np.load('enrichment_baseline_simple_brca.npy')
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[2])



plt.title('Enrichment of GoTerms (x20 rep)')
plt.xlabel('# clusters')
plt.ylabel('enrichment')
plt.legend(['mu10NN', 'mu20NN', 'mu40NN', '3NN', '5NN', '10NN', 'raw spec clus'],loc=4)

#plt.legend(le, ['mu20NN', '5NN', 'raw spec clus', 'treatment', 'control'], loc=2)
plt.ylim([0.9,2.1])




plt.figure(figsize=(7,5))
le = list()

num_clus = np.linspace(5, 100, 20)
treat = np.load('enrichment_baseline_mukNN_'+str(10)+'_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='^', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
treat = np.load('enrichment_baseline_mukNN_'+str(20)+'_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
treat = np.load('enrichment_baseline_mukNN_'+str(40)+'_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
    #plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
#    markersize = 5, color = RGB_tuples[0])

treat = np.load('enrichment_baseline_kNN_'+str(3)+'_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='^', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
treat = np.load('enrichment_baseline_kNN_'+str(5)+'_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
treat = np.load('enrichment_baseline_kNN_'+str(10)+'_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
#plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
#    markersize = 5, color = RGB_tuples[1])

treat = np.load('enrichment_baseline_simple_brca.npy')
treat = [i/k for k,i in enumerate(treat)]
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[2])



plt.title('Enrichment of GoTerms (x20 rep)')
plt.xlabel('# clusters')
plt.ylabel('enrichment')
plt.legend(['mu10NN', 'mu20NN', 'mu40NN', '3NN', '5NN', '10NN', 'raw spec clus'],loc=4)

#plt.legend(le, ['mu20NN', '5NN', 'raw spec clus', 'treatment', 'control'], loc=2)
#plt.ylim([0.9,2.1])
'''
plt.figure(figsize=(7,5))
le = list()

num_clus = np.linspace(5, 100, 20)
treat = np.load('enrichment_baseline_mukNN_'+str(20)+'_brca.npy')
ctrl = np.load('enrichment_baseline_mukNN_'+str(20)+'_ctrl_brca.npy')
ratio[0,:] = 1.0*np.divide(treat,ctrl)
l, = plt.plot([],[], marker='', linestyle='-', \
    color = RGB_tuples[0], linewidth=2)
le.append(l)
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])
plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])

treat = np.load('enrichment_baseline_kNN_'+str(5)+'_brca.npy')
ctrl = np.load('enrichment_baseline_kNN_'+str(5)+'_ctrl_brca.npy')
ratio[1,:] = 1.0*np.divide(treat,ctrl)
l, = plt.plot([],[], marker='', linestyle='-', \
    color = RGB_tuples[1], linewidth=2)
le.append(l)
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])
plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])

treat = np.load('enrichment_baseline_simple_brca.npy')
ctrl = np.load('enrichment_baseline_simple_ctrl_brca.npy')
ratio[2,:] = 1.0*np.divide(treat,ctrl)
l, = plt.plot([],[], marker='', linestyle='-', \
    color = RGB_tuples[2], linewidth=2)
le.append(l)
plt.plot(num_clus, treat, marker='o', linestyle='-', \
    markersize = 5, color = RGB_tuples[2])
plt.plot(num_clus, ctrl, marker='s', linestyle='-', \
    markersize = 5, color = RGB_tuples[2])

l, = plt.plot([],[], marker='o', linestyle='', \
    color = 'k', linewidth=2)
le.append(l)
l, = plt.plot([],[], marker='s', linestyle='', \
    color = 'k', linewidth=2)
le.append(l)

plt.title('Enrichment of GoTerms (x20 rep)')
plt.xlabel('# clusters')
plt.ylabel('enrichment')
plt.legend(le, ['mu20NN', '5NN', 'raw spec clus', 'treatment', 'control'], loc=2)
plt.ylim([0.9,1.7])

'''




'''



plt.figure(figsize=(7,5))

num_clus = np.linspace(5, 100, 20)
plt.plot(num_clus, ratio[0,:], marker='D', linestyle='-', \
    markersize = 5, color = RGB_tuples[0])

l, = plt.plot([],[], marker='D', linestyle='-', \
    color = RGB_tuples[0], linewidth=2)
le.append(l)

plt.plot(num_clus, ratio[1,:], marker='D', linestyle='-', \
    markersize = 5, color = RGB_tuples[1])

l, = plt.plot([],[], marker='D', linestyle='-', \
    color = RGB_tuples[1], linewidth=2)
le.append(l)


plt.plot(num_clus, ratio[2,:], marker='D', linestyle='-', \
    markersize = 5, color = RGB_tuples[2])

l, = plt.plot([],[], marker='D', linestyle='-', \
    color = RGB_tuples[2], linewidth=2)
le.append(l)


plt.title('Enrichment ratio of GoTerms (x20 rep)')
plt.xlabel('# clusters')
plt.ylabel('enrichment ratio: treatment/control')

plt.legend(le, ['mu20NN', '5NN', 'raw spec clus'], loc=2)

'''


'''

plt.figure(figsize=(7,8))

M = 0
csize1 = np.load('classsize_40_baseline_mukNN_'+str(20)+'_brca.npy')
csize2 = np.load('classsize_40_baseline_kNN_'+str(5)+'_brca.npy')
csize3 = np.load('classsize_40_baseline_simple_brca.npy')
M = max(max(csize1), max(csize2), max(csize3))
plt.subplot(3,1,1)
plt.title('mu20NN, 40-means')
plt.hist(csize1,range=(0,M),log=False,bins=M-1, \
    color = RGB_tuples[0])
plt.ylim([0,150])
#plt.xlabel('cluster size')
plt.ylabel('freq')
plt.subplot(3,1,2)
plt.title('5NN, 40-means')
plt.hist(csize2,range=(0,M),log=False,bins=M-1, \
    color = RGB_tuples[1])
plt.ylim([0,150])
#plt.xlabel('cluster size')
plt.ylabel('freq')
plt.subplot(3,1,3)
plt.title('raw graph, 40-means')
plt.hist(csize3,range=(0,M),log=False,bins=M-1, \
    color = RGB_tuples[2])
plt.ylim([0,150])
plt.xlabel('cluster size')
plt.ylabel('freq')

'''

'''

plt.figure(figsize=(7,8))

M = 0
csize1 = np.load('classsize_40_baseline_mukNN_'+str(20)+'_brca.npy')
csize2 = np.load('classsize_40_baseline_kNN_'+str(5)+'_brca.npy')
csize3 = np.load('classsize_40_baseline_simple_brca.npy')
M = max(max(csize1), max(csize2), max(csize3))
plt.subplot(3,1,1)
plt.title('mu20NN, 40-means')
#False
plt.hist(csize1,range=(0,M),log=False,bins=M-1, \
    color = RGB_tuples[0])
plt.ylim([0,150])
plt.xlim([0,30])
#plt.xlabel('cluster size')
plt.ylabel('freq')
plt.subplot(3,1,2)
plt.title('5NN, 40-means')
plt.hist(csize2,range=(0,M),log=False,bins=M-1, \
    color = RGB_tuples[1])
plt.ylim([0,150])
plt.xlim([0,30])
#plt.xlabel('cluster size')
plt.ylabel('freq')
plt.subplot(3,1,3)
plt.title('raw graph, 40-means')
plt.hist(csize3,range=(0,M),log=False,bins=M-1, \
    color = RGB_tuples[2])
plt.ylim([0,150])
plt.xlim([0,30])
plt.xlabel('cluster size')
plt.ylabel('freq')

'''

'''
np.load('classsize_40_baseline_mukNN_'+str(kNN)+'_brca',classsize)

np.load('enrichment_baseline_kNN_'+str(kNN)+'_brca', enrich)
np.load('enrichment_baseline_kNN_'+str(kNN)+'_ctrl_brca', enrich_ctrl)


np.load('classsize_40_baseline_kNN_'+str(kNN)+'_brca',classsize)

np.load('enrichment_baseline_simple_brca', enrich)
np.load('enrichment_baseline_simple_ctrl_brca', enrich_ctrl)

np.load('classsize_40_baseline_simple_brca',classsize)
'''

print 'Done!'