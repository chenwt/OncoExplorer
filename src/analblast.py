# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 10:16:44 2017

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

def readBLAST(filename):
    x,y = [],[]
    f = open(filename, 'r')
    fo = open('out5.txt', 'w')
    for line in f:
        l = line.strip()
        if l.startswith('ENSG'):
            l =l.split()
            
            sc = float(l[1])
#            if 'e' in sc:
#                sc = sc.split('e')
#                sc0, sc1 = float(sc[0]), int(sc[1])
#                sc = 1.0*sc0**sc1
#            else:
#                sc = float(sc)
            ev = float(l[2])
            if ev == 0: ev = 1e-180
            #print >> fo, l
            print >> fo, sc, ev
            x.append(sc)
            y.append(ev)
        #print l
    
    fo.close()
    f.close()
    x,y = np.asarray(x), np.asarray(y)
    xl = np.log10(x)
    yl = np.log10(y)
    return xl, yl


plt.figure(figsize=(10,5))

path = '../blastresult/out_unmask.txt'
xl, yl = readBLAST(path)
plt.plot(xl,yl,linestyle='',marker='o',markersize=0.5,color='k')
plt.xlabel('log10(Bit-score)')
plt.ylabel('log10(E-value)')

#path = '../blastresult/out_hard_mask_mask.txt'
#xl, yl = readBLAST(path)
#plt.plot(xl,yl,linestyle='',marker='o',markersize=2,color='b')
#
#path = '../blastresult/out_hard_dust_mask.txt'
#xl, yl = readBLAST(path)
#plt.plot(xl,yl,linestyle='',marker='o',markersize=2,color='g')






#plt.plot(2,yl[1],linestyle='',marker='o',markersize=5)

