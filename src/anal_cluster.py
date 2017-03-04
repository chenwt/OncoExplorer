# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 10:04:45 2017

@author: Yifeng Tao
"""

from collections import defaultdict as dd
import numpy as np
import scipy.io
import os
from random import shuffle
import math


set_g = set()
sgaAsga2pp = dd(float)
path = 'sgaAsga2pp.txt'
f = open(path, 'r')
for line in f:
    l = line.strip().split('\t')
    sga1, sga2, prob = l[0], l[1], float(l[2])
    sgaAsga2pp[(sga1, sga2)] = prob
    set_g.add(sga1)
    set_g.add(sga2)

f.close()


#count = 0
y2g = dd(set)
path = 'cluster.txt'
f = open(path, 'r')
for line in f:
    l = line.strip().split('\t')
    g, y = l[0], int(l[1])
    #count += 1
    if g in set_g:
        
        #print count
        y2g[y].add(g)
f.close()

go2prob = dd(float)
path = 'go2prob.txt'
f = open(path, 'r')
for line in f:
    l = line.strip().split('\t')
    go, prob = l[0], float(l[1])
    go2prob[go] = prob
f.close()

g2go = dd(set)
path = 'g2go.txt'
f = open(path, 'r')
for line in f:
    l = line.strip().split('\t')
    g, go = l[0], l[1]
    g2go[g].add(go)
f.close()

for y in y2g.keys():
    for it, g in enumerate(y2g[y]):
        if it == 0:
            u = g2go[g]
        else:
            u = u.intersection(g2go[g])
    mP = 1.0
    mGo = ''
    for go in u:
        v = go2prob[go]
        if v <= mP:
            mP = v
            mGo = go
    print 'class: '+str(y)+'\tsize='+str(len(y2g[y]))+'\tminGo='+mGo+'\tminP(go)='+str(mP)
    #print len(u)
    

path = 'c2gene2.txt'
f = open(path, 'w')
for y in y2g.keys():
    for g in y2g[y]:
        print >> f, str(y)+'\t'+g
f.close()

path = 'c2gene3.txt'
f = open(path, 'w')
for y in y2g.keys():
    print >> f, str(y)
    for g in y2g[y]:
        print >> f, g
f.close()







