# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 04:52:07 2017

@author: Yifeng Tao
"""

# Generate mutual kNN graph.

import numpy as np
import numpy.linalg as linalg
from collections import defaultdict as dd
import matplotlib.pyplot as plt
from random import shuffle
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.manifold import TSNE


itertime = 10
KM = np.array([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
enrich = np.zeros([len(KM),],dtype = float)

for iter in range(itertime):
    
    
    kNN = 40
    sga2sgaAaff = dd(set)
    sgaAsga2aff = dd(int)
    set_gene = set()
    path = 'sga2sgaAff_baseline_brca.txt'
    f = open(path, 'r')
    pathout = 'out.txt'
    fo = open(pathout, 'w')
    next(f)
    for line in f:
        l = line.strip().split('\t')
        sga1, sga2, aff = l[0], l[1], int(l[2])
        sga2sgaAaff[sga1].add((sga2, aff))
        set_gene.add(sga1)
        set_gene.add(sga2)
    f.close()
    
    set_gene = list(set_gene)
    
    for sga in set_gene:
        context = sga2sgaAaff[sga]
        context = list(context)    
        # shuffle
        
        for i in range(kNN):
            shuffle(context)
        
            Maff = 0.0
            Mgene = 'helloword'
            for line in context:
                sgacontext, aff = line[0], line[1]
                if aff > Maff:
                    Maff = aff
                    Mgene = sgacontext
            sgaAsga2aff[(sga, Mgene)] = aff
            #sgaAsga2aff[(Mgene, sga)] = aff
            #print sga, Mgene, Maff
            context.remove((Mgene, Maff))
            if len(context) == 0: break
        
        #sgaAsga2aff[(sga1, sga2)] = aff
    
        #print >> fo, l
    
    gene2id = dd()
    id2gene = dd()
    
    k = 0
    for g in set_gene:
        gene2id[g] = k
        id2gene[k] = g
        k += 1
    
    
    len_gene = len(set_gene)
    W = np.zeros([len_gene, len_gene], float)
    for i in range(len_gene):
        for j in range(len_gene):
            W[i][j] = sgaAsga2aff[(id2gene[i], id2gene[j])]*sgaAsga2aff[(id2gene[j], id2gene[i])]
    
    
    #D = W.sum(axis=0)
    #Dinv = np.diag(1.0/D)
    #
    #D = np.diag(D);
    #L =  D - W
    #
    #lambda0, V0 = linalg.eig(L)
    #
    #ordr =lambda0.argsort()
    #lam = lambda0[ordr]
    #V = V0[:,ordr]
    #
    #
    #
    #plt.plot(lam)
    #plt.ylabel('labmda')
    #plt.xlabel('k')
    #plt.title('SGA in BRCA')
    #plt.ylim((0,50))
    
    
    
    
    
    path = 'sga2sgaAff_baseline_mukNN_'+str(kNN)+'_brca.txt'
    f = open(path, 'w')
    print >> f, 'Source\tTarget\tWeight'
    for line in sgaAsga2aff.keys():
            sga1, sga2 = line[0], line[1]
            overlap = sgaAsga2aff[(sga1, sga2)]*sgaAsga2aff[(sga2, sga1)]
            #sga2, overlap = line[0], line[1]
            if overlap == 0: continue
            if sga1 == sga2: continue
            print >> f, sga1 + '\t' + sga2 + '\t' + str(overlap)
            #print >> f, sga1 + '\t' + sga2 + '\t' + overlap
    f.close()
    
    
    
    
    
    
    sga2sgaAaff = dd(set)
    set_gene = set()
    
    path = 'sga2sgaAff_baseline_mukNN_'+str(kNN)+'_brca.txt'
    f = open(path, 'r')
    next(f)
    for line in f:
        l = line.strip().split('\t')
        sga1, sga2, aff = l[0], l[1], int(l[2])
        sga2sgaAaff[sga1].add((sga2, aff))
        set_gene.add(sga1)
        set_gene.add(sga2)
    f.close()
    
    set_gene = list(set_gene)
    
    gene2id = dd()
    id2gene = dd()
    
    k = 0
    for g in set_gene:
        gene2id[g] = k
        id2gene[k] = g
        k += 1
    
    
    len_gene = len(set_gene)
    W = np.zeros([len_gene, len_gene], float)
    for i in range(len_gene):
        for j in range(len_gene):
            W[i][j] = sgaAsga2aff[(id2gene[i], id2gene[j])]*sgaAsga2aff[(id2gene[j], id2gene[i])]
    
    D = W.sum(axis=0)
    Dinv = np.diag(1.0/D)
    
    D = np.diag(D);
    L =  D - W
    #L = D - W;
    L = np.dot(Dinv,L)
    
    #print 2
    
    lambda0, V0 = linalg.eig(L)
    
    #print 3
    #L = diag(1./sum(W))*L;
    #[V0, lambda0] = eig(L);
    ordr =lambda0.argsort()
    lam = lambda0[ordr]
    #V = V0
    V = V0[:,ordr]

    
    
    gene2goid = dd(set)
    path = 'gene2goId.txt'
    f = open(path, 'r')
    #next(f)
    for line in f:
        l = line.strip().split('\t')
        g, goid = l[0], l[1]
        gene2goid[g].add(goid)
    f.close()
    
    
    
    #KM = np.array([40,50])
    
    t = 0
    for km in KM:
        #km = 50
        X = V[:,0:km]
        
        #kmeans = KMeans(n_clusters=km, random_state=0).fit(X)
        kmeans = KMeans(n_clusters=km).fit(X)
        y = kmeans.labels_
        
        
        
        len_gene = len(set_gene)
        
        total_u = 0.0
        total_d = 0.0
        within_u = 0.0
        within_d = 0.0
        for i in range(len_gene):
            g1 = id2gene[i]
            for j in range(len_gene):
                g2 = id2gene[j]
                l1 = y[i]
                l2 = y[j]
                total_d += 1
                if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                    total_u += 1
                if y[i] == y[j]:
                    within_d += 1
                    if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                        within_u += 1
        #print within_d, total_d, within_u, total_u
        print km, 1.0*within_u/within_d*total_d/total_u
        enrich[t] += 1.0*within_u/within_d*total_d/total_u
        t += 1
    

    

t = 0
km = 40
#km = 50
X = V[:,0:km]

kmeans = KMeans(n_clusters=km, random_state=0).fit(X)
y = kmeans.labels_



len_gene = len(set_gene)

total_u = 0.0
total_d = 0.0
within_u = 0.0
within_d = 0.0
for i in range(len_gene):
    g1 = id2gene[i]
    for j in range(len_gene):
        g2 = id2gene[j]
        l1 = y[i]
        l2 = y[j]
        total_d += 1
        if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
            total_u += 1
        if y[i] == y[j]:
            within_d += 1
            if len(gene2goid[g1].intersection(gene2goid[g2])) > 0:
                within_u += 1
#print within_d, total_d, within_u, total_u
print km, 1.0*within_u/within_d*total_d/total_u


#X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
#y = np.array([0,0,0,1,1,1])
#plt.figure()
#h = X[:,0]
#v = X[:,1]
#plt.plot(h,v, 'ro')



pca = PCA(n_components=2)
Xnew = pca.fit_transform(X)
h = Xnew[:,0]
v = Xnew[:,1]


#plt.figure()
#plt.plot(h,v, 'ro')

df = pd.DataFrame(dict(x=h, y=v, label=y))

groups = df.groupby('label')

plt.figure(figsize=(8, 8))
#fig, ax = plt.subplots()
for name, group in groups:
    plt.plot(group.x, group.y, marker='o', linestyle='', markersize = 5, label=name)
#ax.legend()
plt.xlim(-1,1)
plt.ylim(-1,1)


#X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])


model = TSNE(n_components=2, random_state=0)
#np.set_printoptions(suppress=True)
Xnew = model.fit_transform(X)
h = Xnew[:,0]
v = Xnew[:,1]


#plt.figure()
#plt.plot(h,v, 'ro')

df = pd.DataFrame(dict(x=h, y=v, label=y))

groups = df.groupby('label')

#ax = plt.figure(figsize=(8, 8))
fig, ax = plt.subplots(figsize=(8, 8))




#set_gene = list(set_gene)



for i in range(len(set_gene)):
    sga1 = set_gene[i]
    for line in sga2sgaAaff[sga1]:
        sga2, aff = line[0], line[1]
        sga1id = gene2id[sga1]
        sga2id = gene2id[sga2]
        ax.plot([h[sga1id], h[sga2id]], [v[sga1id], v[sga2id]], marker='.', linestyle='-', color = [0.9, 0.9, 0.9] )

for name, group in groups:
    ax.plot(group.x, group.y, marker='o', linestyle='', markersize = 5, label=name)
#ax.legend()
    
    
    
threshold = 300

blacklist = set()
path = 'degree_SGA_SGA_.txt'
f = open(path, 'r')
next(f)
for line in f:
    l = line.strip().split('\t')
    sga, degree = l[0], int(l[1])
    if degree > threshold: blacklist.add(sga)
f.close()
#print 'len(blacklist)=', len(blacklist)
#print blacklist
for sga in blacklist:
    if sga not in gene2id.keys(): continue
    sgaid = gene2id[sga]
    ax.annotate(sga, xy=(h[sgaid], v[sgaid])#, xytext=(3, 1.5),
            #arrowprops=dict(shrink=0.05),
            )
    ax.plot(h[sgaid], v[sgaid], marker='*')


plt.title('t-SNE of spectral clustering (k = '+str(km)+')')
#plt.xlim(-1,1)
#plt.ylim(-1,1)

plt.figure()
plt.plot(lam)
plt.ylabel('labmda')
plt.xlabel('k')
plt.title('SGA in BRCA')
plt.ylim((0,2))


plt.figure()
plt.plot(KM, enrich/itertime)
plt.xlabel('k')
plt.ylabel('enrichment')

np.save('enrichment_baseline_mukNN_'+str(kNN)+'_brca', enrich/itertime)
plt.title('averaged over '+str(itertime)+' replications')


print 'Done!'













