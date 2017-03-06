# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 06:51:38 2017

@author: Yifeng Tao
"""

# Compare pairwise distance.
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

choice = 14
if choice == 1:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'Baseline1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.5*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    ns = len(aff_goterm1)
    nsx = 0.00*np.random.normal(size=(ns,))
    nsy = 0.00*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(Baseline affinity)')
    plt.xlabel('GoTerm affinity')
    #plt.title('SGA in BRCA')
    #plt.ylim((0,2))
if choice == 2:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'Baseline1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.5*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(Baseline affinity)')
    plt.xlabel('log10(GoTerm affinity)')
    #plt.title('SGA in BRCA')
    #plt.ylim((0,2))    

if choice == 3:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'Baseline1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'WeightedOverlap_max.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.5*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(Baseline affinity)')
    plt.xlabel('log10(Weighted overlap)')
    #plt.title('SGA in BRCA')
    #plt.ylim((0,2))    


if choice == 4:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'WeightedOverlap_max.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.5*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.00*np.random.normal(size=(ns,))
    nsy = 0.00*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(weighted overlap)')
    plt.xlabel('log10(GoTerm affinity)')
    #plt.title('SGA in BRCA')
    #plt.ylim((0,2))    

if choice == 5:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'LuAff.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.00*np.random.normal(size=(ns,))
    nsy = 0.00*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(Lu affinity)')
    plt.xlabel('log10(GoTerm affinity)')


if choice == 6:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'TFIDF_max.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.005*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.00*np.random.normal(size=(ns,))
    nsy = 0.00*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(TF-IDF)')
    plt.xlabel('log10(GoTerm affinity)')


if choice == 7:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'Jaccard.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('log10(Jaccard)')
    plt.xlabel('log10(GoTerm affinity)')

if choice == 8:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'spec.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    #aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    #aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('Spectral embedding')
    plt.xlabel('log10(GoTerm affinity)')

if choice == 9:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'spec.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'Baseline1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    #aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    #aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.00*np.random.normal(size=(ns,))
    nsy = 0.00*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('Spectral embedding')
    plt.xlabel('log10(Baseline)')


if choice == 10:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'specmukNN.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    #aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    #aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('Spectral embedding(mukNN)')
    plt.xlabel('log10(GoTerm affinity)')

if choice == 11:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'specmukNN.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'Baseline1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    #aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    #aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.00*np.random.normal(size=(ns,))
    nsy = 0.00*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('Spectral embedding(mukNN)')
    plt.xlabel('log10(Baseline)')


if choice == 12:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'LuDist.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    aff_baseline1 = [0.0001*(line==0)+line*(line!=0) for line in aff_baseline1]
    aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    #plt.ylabel('Spectral embedding(mukNN)')
    #plt.xlabel('log10(Baseline)')
    plt.ylabel('log10(Lu dist)')
    plt.xlabel('log10(GoTerm affinity)')
    #BUGs in zeros level.

if choice == 13:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'DWsga_wl2_nw200.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'Baseline1.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    #aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    #aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('DWsga')
    plt.xlabel('log10(Baseline)')

if choice == 14:
    sgaAsga2aff_baseline1 = dd(float)
    path = 'pairdist/'+'brca_'+'DW_wl2_nw200.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_baseline1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    sgaAsga2aff_goterm1 = dd(float)
    path = 'pairdist/'+'brca_'+'goterm.txt'
    f = open(path, 'r')
    for line in f:
        l = line.strip().split('\t')
        sgaAsga2aff_goterm1[(l[0],l[1])] = float(l[2])
        
    f.close()
    
    print 'calculating...'
    set_line = set(sgaAsga2aff_baseline1.keys()).intersection(set(sgaAsga2aff_goterm1.keys()))
    set_line = list(set_line)
    aff_goterm1 = [sgaAsga2aff_goterm1[line] for line in set_line]
    aff_baseline1 = [sgaAsga2aff_baseline1[line] for line in set_line]
    #aff_baseline1 = [0.0003*(line==0)+line*(line!=0) for line in aff_baseline1]
    #aff_baseline1 = [np.log10(line) for line in aff_baseline1]
    aff_goterm1 = [0.5*(line==0)+line*(line!=0) for line in aff_goterm1]
    aff_goterm1 = [np.log10(line) for line in aff_goterm1]
        
    ns = len(aff_goterm1)
    nsx = 0.01*np.random.normal(size=(ns,))
    nsy = 0.01*np.random.normal(size=(ns,))
    
    print 'plotting...'
    plt.figure(figsize=(5,5))
    plt.plot(aff_goterm1+nsx,aff_baseline1+nsy,marker='.', linestyle='', color = 'k' )
    #plt.xlim([-0.5,1.5])
    plt.grid()
    plt.ylabel('DW')
    plt.xlabel('log10(goterm)')




