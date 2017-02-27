# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 15:49:15 2017

@author: Yifeng Tao
"""

# myDijkstra.py

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

def getGGaff(cancer):
    #set_g = set()
    sgaAsga2aff = dd(float)    
    sga2deg = dd(set)

    flag = True
    f = open('../ensemble.txt', 'r')
    next(f)
    k = 0
    for line in f:
        line = line.strip().split('\t')
        patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
        #print can, sga, deg, prob
        #if sga not in set_g: continue
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
    set_g = sga2deg.keys()
    for sga1 in sga2deg.keys():
        #traversedset.add(sga1)
        
        for sga2 in sga2deg.keys(): 
            if sga1 == sga2: continue
            #print sga1, sga2
            overlap = 1.0*len(sga2deg[sga1].intersection(sga2deg[sga2]))
            if overlap == 0: overlap = 0.5
            sgaAsga2aff[(sga1,sga2)] = overlap
    return sgaAsga2aff, set_g




#print 'Read our dataset...'
#sgaAsga2aff, set_g = getGGaff('pancan')



print 'Read into dawnrank...'
distances = dd(set)

set_G = set()
path = 'graph_dawnrank.txt'
f = open(path, 'r')
for t, line in enumerate(f):
    l = line.strip().split('\t')
    if l[0] == 'resultFrom': continue
    src = l[1].lower()
    dst = l[2].lower()
    set_G.add(src)
    set_G.add(dst)
    #if distance[src] == 0:
    #    distance[src] = dd()
    #if distance[dst] == 0:
    #    distance[dst] = dd()
    distances[src].add(dst)
    distances[dst].add(src)
    #print t, src, dst
    #if t >= 100: break
    
f.close()





print 'Dijkstra...'


#set_g = set(set_g)

nodes = set_G

nodes = list(nodes)
#nodes = ('A', 'B', 'C', 'D')

#nodes =  set_g


#nodes = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'V')
#nodes = ('A', 'B', 'C', 'H', 'V')
#distances = {
#    'B': {'A': 1, 'D': 1, 'G': 1},
#    'A': {'B': 1, 'D': 1, 'E': 1, 'F' :1},
#    'D': {'B': 1, 'G': 1, 'E': 1, 'A': 1},
#    'G': {'B': 1, 'D': 1, 'C': 1},
#    'C': {'G': 1, 'E': 1, 'F': 1},
#    'E': {'A': 1, 'D': 1, 'C': 1, 'F': 1},
#    'F': {'A': 1, 'E': 1, 'C': 1}}
#distances = {
#    'B': {'A', 'D', 'G'},
#    'A': {'B', 'D', 'E', 'F'},
#    'D': {'B', 'G', 'E', 'A'},
#    'G': {'B', 'D', 'C'},
#    'C': {'G', 'E', 'F'},
#    'E': {'A', 'D', 'C', 'F'},
#    'F': {'A', 'E', 'C'},
#    'H': {'V'},
#    'V': {'H'}}
#distances = {
#    'A': {'B'},
#    'B': {'A', 'C'},
#    'C': {'B'},
#    'H': {'V'},
#    'V': {'H'}}   
#nodes = list(set_G)

path = 'shortest.txt'
f = open(path, 'w')

count = 0
for current in nodes[1:3]:
    
    if count % 5 == 0:
        print count, current
    count += 1
    s = current
    unvisited = {node: 999999999999 for node in nodes} #using None as +inf
    #break
    visited = {}
    #current = 'B'
    currentDistance = 0
    unvisited[current] = currentDistance
    
    while True:
        for neighbour in distances[current]:
            distance = 1
        #for neighbour, distance in distances[current].items():
            if neighbour not in unvisited: continue
            newDistance = currentDistance + distance
            if unvisited[neighbour] == 999999999999 or unvisited[neighbour] > newDistance:
                unvisited[neighbour] = newDistance
        visited[current] = currentDistance
        del unvisited[current]
        if not unvisited: break
        #candidates = [node for node in unvisited.items() if node[1]]
        candidates = [node for node in unvisited.items()]
        current, currentDistance = sorted(candidates, key = lambda x: x[1])[0]
    
    for sga2 in visited:
        print >> f, s+'\t'+sga2+'\t'+str(visited[sga2])
    #print(visited)

f.close()

#EOF.