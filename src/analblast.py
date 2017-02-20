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

def getGGaff(cancer, set_g):

    sgaAsga2aff = dd(float)    
    sga2deg = dd(set)

    flag = True
    f = open('../src/ensemble.txt', 'r')
    next(f)
    k = 0
    for line in f:
        line = line.strip().split('\t')
        patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
        #print can, sga, deg, prob
        if sga not in set_g: continue
        if cancer == 'pancan': # pancan is alway true
            flag = True
        elif can == cancer: # cancer name match specific can
            flag = True
        else: # cancer does not match the specific can
            flag = False
        if not flag: continue
        k = k+1
        
        sga2deg[sga].add(deg)

        #print len(sga2deg.keys())
        
    for sga1 in sga2deg.keys():
        #traversedset.add(sga1)
        for sga2 in sga2deg.keys(): 
            if sga1 == sga2: continue
            #print sga1, sga2
            overlap = 1.0*len(sga2deg[sga1].intersection(sga2deg[sga2]))
            if overlap == 0: overlap = 0.5
            sgaAsga2aff[(sga1,sga2)] = overlap
    return sgaAsga2aff


def readBLAST(filename, e2g):
    sgaAsga2aff = dd(float)
    #x,y = [],[]
    f = open(filename, 'r')
    fo = open('out5.txt', 'w')
    #flag = False
    for line in f:
        l = line.strip()
        if l.startswith('Query='):
            sga1 = l.split('|')[1]
            #print >> fo, sga1
#            flag = True
#        elif l.startswith('Effective'):
#            flag = False
        elif l.startswith('ENSG'):
            l =l.split()
            sga2 = l[0]
            if sga1 == sga2: continue
            sc = float(l[1])
            ev = float(l[2])
            if ev == 0: ev = 1e-180
            #print >> fo, l
            print >> fo, sga1, sga2, sc, ev
            scl, evl = np.log10(sc), np.log10(ev)
            sgaAsga2aff[(e2g[sga1],e2g[sga2])] =  scl
            sgaAsga2aff[(e2g[sga2],e2g[sga1])] =  scl
            #x.append(sc)
            #y.append(ev)
        #print l
    
    fo.close()
    f.close()
#    x,y = np.asarray(x), np.asarray(y)
#    xl = np.log10(x)
#    yl = np.log10(y)
#    
    return sgaAsga2aff




for cancer in ['pancan', 'brca', 'gbm', 'ov']:
    for mymode in ['out_unmask', 'out_hard_dust_mask', 'out_hard_mask_mask']:

#mymode = 'out_unmask'

#cancer = 'brca'
#cancer = 'pancan'
#cancer = 'gbm'
#cancer = 'ov'



        g2e = dd()
        e2g = dd()
        path = '../src/gene2ensembl.txt'
        f = open(path, 'r')
        next(f)
        for line in f:
            l = line.strip().split('\t')
            g, e = l[0], l[1]
            g2e[g] = e
            e2g[e] = g
            #set_g.add(g)
            #print g,e
        f.close()
        
        
        # set of genes that are searched.
        set_g = set()
        
        path = '../src/set_ensembl.txt'
        f = open(path, 'r')
        for line in f:
            l = line.strip()
            set_g.add(e2g[l])
        f.close()
        
        
        
        
        print 'Reading BLAST...'
        
        path = '../blastresult/'+mymode+'.txt'
        sgaAsga2affB = readBLAST(path, e2g)
        
        
        print 'Reading ensemble...'
        
        sgaAsga2aff = getGGaff(cancer, set_g)
        
        
        #set_ggB = sgaAsga2affB.keys()
        
        print 'Calculating...'
        
        set_gg = sgaAsga2aff.keys()
        set_gg = list(set_gg)
        xd = [sgaAsga2affB[gg] for gg in set_gg]
        yd = [sgaAsga2aff[gg] for gg in set_gg]
        xd = np.asarray(xd)
        yd = np.asarray(yd)
        
        print 'Saving...'
        np.save('xd_'+cancer+'_'+mymode, xd)
        np.save('yd_'+cancer+'_'+mymode, yd)
        
        '''
        #print 'Ploting...'
        plt.figure(figsize=(5,5))
        #for gg in set_gg:
        plt.plot(xd,yd,linestyle='',marker='o',markersize=2,color='b')
        
        plt.xlabel('log10(BLAST Bit-score)')
        plt.ylabel('Baseline affinity')
        
        #plt.xlabel('log10(Bit-score)')
        #plt.ylabel('log10(E-value)')
        '''
        
        
        '''
        xl, yl, sgaAsga2affB = readBLAST(path)
        plt.figure(figsize=(5,5))
        plt.plot(xl,yl,linestyle='',marker='o',markersize=0.5,color='k')
        plt.xlabel('log10(Bit-score)')
        plt.ylabel('log10(E-value)')
        
        
        '''
        #path = '../blastresult/out_hard_mask_mask.txt'
        #xl, yl = readBLAST(path)
        #plt.plot(xl,yl,linestyle='',marker='o',markersize=2,color='b')
        #
        #path = '../blastresult/out_hard_dust_mask.txt'
        #xl, yl = readBLAST(path)
        #plt.plot(xl,yl,linestyle='',marker='o',markersize=2,color='g')
        
        
        
        
        
        
        #plt.plot(2,yl[1],linestyle='',marker='o',markersize=5)



#EOF.