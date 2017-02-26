# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 10:56:00 2017

@author: Yifeng Tao
"""

# Prepare piptide data

path = 'raw.fasta'
f = open(path, 'r')


i = 0
#flag = 'seq'
sequence = ''
seq = ''
for line in f:
    l = line.strip()
    if l.startswith('>ENSG'):
        if seq != '':
            i += 1
        #i += 1
        
        #flag = 'ensbl'
        sequence += seq
        seq = l.split('|')[0]+'i'+str(i)+'\n'
    else:
        if l.startswith('Sequence unavailable'):
            seq = ''
            continue
        seq += l.strip('*')+'\n'
        
    #print l
f.close()
patho = 'pptd.fa'
fo = open(patho, 'w')
print >> fo, sequence

fo.close()















#EOF.