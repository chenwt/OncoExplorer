# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 05:04:23 2017

@author: Yifeng Tao
"""

# Examine the result of fasta
path = 'seq.fasta'
f = open(path, 'r')
pathout = 'myseq.fa'
fo = open(pathout, 'w')
for line in f:
    l = line.strip()
    if l.startswith('>'):
        l = l.split('|')[0][1:]
        print >> fo, '>lcl|'+l
    else:
        print >> fo, l

f.close()
fo.close()

print 'Done!'