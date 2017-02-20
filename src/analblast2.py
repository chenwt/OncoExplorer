# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 04:07:21 2017

@author: Yifeng Tao
"""
# Analyze the result from BLAST.

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

for cancer in ['pancan', 'brca', 'gbm', 'ov']:
    for mymode in ['out_unmask', 'out_hard_dust_mask', 'out_hard_mask_mask']:
#for cancer in ['pancan']:
#    for mymode in ['out_unmask']:
        xd = np.load('xd_'+cancer+'_'+mymode+'.npy')
        xd = [1.0*(i==0)+i for i in xd]
        yd = np.load('yd_'+cancer+'_'+mymode+'.npy')
        yd = np.log10(yd)
        plt.figure(figsize=(5,5))
        #for gg in set_gg:
        plt.plot(xd,yd,linestyle='',marker='o',markersize=1,color='k')
        
        plt.xlabel('log10(BLAST Bit-score)')
        plt.ylabel('log10(Baseline affinity)')
        plt.xlim([0,5])
        plt.ylim([-0.5,3])
        plt.title(cancer+':BLAST '+mymode[4:])
        