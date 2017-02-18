# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 04:13:47 2017

@author: Yifeng Tao
"""

from collections import defaultdict as dd

# Prepare data for deepwalk.

set_sga = set()
set_deg = set()
sga2deg = dd(set)
deg2sga = dd(set)


path = '../src/ensemble.txt'
f = open(path, 'r')
next(f)
k = 0
for line in f:
    l = line.strip().split('\t')
    can, sga, deg = l[2], l[4], l[5]
    deg = deg+'_deg'
    if can == 'brca':
        set_sga.add(sga)
        set_deg.add(deg)
        sga2deg[sga].add(deg)
        deg2sga[deg].add(sga)
        #print can, sga, deg
    k += 1
    #if k == 10: break
f.close()

remain_deg = set()
for deg in set_deg:
    remain_deg.add(deg)
    
for deg in set_deg:
    if len(deg2sga[deg]) == 1:
        remain_deg.remove(deg)

for sga in set_sga:
    sga2deg[sga] = sga2deg[sga].intersection(remain_deg)

set_deg = remain_deg


set_sga = list(set_sga)
set_deg = list(set_deg)
gene2id = dd(int)
k = 1
for sga in set_sga:
    gene2id[sga] = k
    k += 1
for deg in set_deg:
    gene2id[deg] = k
    k += 1

f = open('LUT_DW_brca.txt', 'w')
for sga in set_sga:
    print >> f, sga+'\t'+str(gene2id[sga])
for deg in set_deg:
    print >> f, deg+'\t'+str(gene2id[deg])
f.close()

f = open('graph_brca.adjlist', 'w')
for sga in set_sga:
    s = str(gene2id[sga])
    for deg in sga2deg[sga]:
        s = s + ' ' + str(gene2id[deg])
    print >> f, s
for deg in set_deg:
    s = str(gene2id[deg])
    for sga in deg2sga[deg]:
        s = s + ' ' + str(gene2id[sga])
    print >> f, s
f.close()
















#EOF.