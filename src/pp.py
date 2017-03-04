# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 06:24:49 2017

@author: Yifeng Tao
"""

#Get two files:
#goTerm frequency of genes.

from collections import defaultdict as dd




set_g = set()
path = 'degree_SGA_SGA_.txt'
f = open(path, 'r')
next(f)
for line in f:
    l = line.strip().split('\t')[0]
    set_g.add(l)
    #print l
f.close()

graph = dd(set)

path = 'goAnn.cfacts'
f = open(path, 'r')
for line in f:
    l = line.strip().split('\t')
    g = l[1]
    go = l[2]
    graph[g].add(go)
    
f.close()

path = 'isA.cfacts'
f = open(path, 'r')
for line in f:
    l = line.strip().split()
    src = l[1]
    dst = l[2]
    graph[src].add(dst)
f.close()


def searchGo(g, graph):
    # BFS
    go = set()
    st = list()
    st.append(g)
    while True:
        if len(st) == 0: break
        curnt = st.pop(0)
        for tp in graph[curnt]:
            st.append(tp)
            go.add(tp)
    return go


g2go = dd(set)
set_g = list(set_g)
l_g = len(set_g)
for g in set_g:
    go = searchGo(g, graph)
    g2go[g] = go

go2prob = dd(float)
path = 'g2go.txt'
l_r = 0
set_r = set()
f = open(path, 'w')
for g in set_g:
    if len(g2go[g]) >= 1:
        #l_r += 1
        set_r.add(g)
    for go in g2go[g]:
        print >> f, g+'\t'+go
        go2prob[go] += 1.0
f.close()
l_r = len(set_r)

path = 'go2prob.txt'
f = open(path, 'w')
for go in go2prob.keys():
    go2prob[go] = 1.0*go2prob[go]/l_r
    print >> f, go+'\t'+str(go2prob[go])
f.close()



sgaAsga2pp = dd(float)
for sga1 in set_r:
    for sga2 in set_r:
        if sga1 == sga2: continue
        sh_go = g2go[sga1].intersection(g2go[sga2])
        if len(sh_go) == 0: continue
        mpp = 1.0        
        for go in sh_go:
            tmp = go2prob[go]
            tmp = 1.0*tmp*tmp
            if tmp < mpp:
                mpp = tmp
        sgaAsga2pp[(sga1,sga2)] = mpp
        sgaAsga2pp[(sga2,sga1)] = mpp

path = 'sgaAsga2pp.txt'
f = open(path, 'w')
for line in sgaAsga2pp.keys():
    sga1, sga2 = line[0], line[1]
    pp = sgaAsga2pp[(sga1, sga2)]
    print >> f, sga1+'\t'+sga2+'\t'+str(pp)

f.close()

path = 'genewithGO.txt'
f = open(path, 'w')
for sga in set_r:
    print >> f, sga
f.close()
















#EOF.