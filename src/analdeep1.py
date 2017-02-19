# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 00:41:25 2017

@author: Yifeng Tao
"""

# analdeep1.py
# plot figure for analdeep.py
#import matplotlib.patches as mpatches
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


data=np.load('true_km20_brca.npy')
data_f=np.load('fake_km20_brca.npy')
NW = [10,20,50,100,200]
WL = [1,2,3,5,10,20,40]
WL = ['walk_len='+str(wl) for wl in WL]
fig = plt.figure(figsize=(8, 8))


RGB_tuples = _get_colors(len(WL))

le = []
for count_wl,wl in enumerate(WL):
    l, = plt.plot([],[], linestyle='-', linewidth=2,color=RGB_tuples[count_wl], label=WL[count_wl])
    plt.plot(NW, data[count_wl,:], linestyle='-', linewidth=2.0, marker='o', markersize=10, color=RGB_tuples[count_wl])
    le.append(l)
    #print count_wl
for count_wl,wl in enumerate(WL):
    plt.plot(NW, data_f[count_wl,:], linestyle='-', linewidth=2, marker='s', markersize=10, color=RGB_tuples[count_wl])

l, = plt.plot([],[],linestyle='-',color='k', marker='o')
le.append(l)
WL.append('treatment')

l, = plt.plot([],[],linestyle='-',color='k', marker='s')
le.append(l)
WL.append('control')
    #plt
plt.title('Enrichment of GoTerm (Deep walk), 20-means x10')
plt.xlabel('#walks/#nodes')
plt.ylabel('enrichment')
plt.xlim([0.0,210])
plt.ylim([1.0,2.0])
plt.legend(le, WL,loc=1)



#line_up, = plt.plot([1,2,3], label='Line 2')
#line_down, = plt.plot([3,2,1], label='Line 1')
#plt.legend([line_up, line_down], ['Line Up', 'Line Down'])