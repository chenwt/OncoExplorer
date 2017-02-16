# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 03:21:50 2017

@author: Yifeng Tao
"""

# LUT of gene name and ensembl id

import mygene
from collections import defaultdict as dd


    #print line


set_gene = set()

path = 'ensemble.txt'
f = open(path, 'r')
next(f)
for line in f:
    l = line.strip().split('\t')
    sga = l[4]
    #print sga
    set_gene.add(sga)
    
f.close()

#print len(set_gene)
set_gene = list(set_gene)


mg = mygene.MyGeneInfo()
#xli = ['tp53', 'brca1']
xli = set_gene
print len(xli)
out = mg.querymany(xli, scopes='symbol', fields='ensembl.gene', species='human')

gene2enid = dd(set)
for line in out:
    #print line
    gene = line[u'query']
    if u'ensembl' not in line.keys():
        #print 
        #print >> f, gene + ' No ensembl!'
        continue
    enid = line[u'ensembl']
    if len(enid) > 1:
        #print >> f, 'Not consider'
        continue
    enid = enid[u'gene']
    # 498
    gene2enid[gene].add(enid)
    #print >> f, gene, enid

path = 'gene2ensembl.txt'
f = open(path, 'w')
path2 = 'set_ensembl.txt'
f2 = open(path2, 'w')
print >> f, 'Gene\tEnsemblId'
for g in gene2enid.keys():
    enid = gene2enid[g]    
    if len(enid) > 1:
        continue
    for enid1 in enid:
        print >> f, g+'\t'+enid1
        print >> f2, enid1
        #print g, enid
 
f.close()
f2.close() 
print 'Done'
#EOF.