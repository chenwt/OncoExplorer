# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd
import numpy as np
import scipy.io
import os
from random import shuffle

def getAffinityBaseline1(cancer):
    sga2sgaAaff = dd(set)
    sga2deg = dd(set)
    k = 0
    flag = True
    with open('/usr1/public/yifeng/Github/inputData/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
            if cancer == 'pancan': # pancan is alway true
                flag = True
            elif can == cancer: # cancer name match specific can
                flag = True
            else: # cancer does not match the specific can
                flag = False
            if not flag: continue
            k = k+1

            sga2deg[sga].add(deg)

    f = open('/usr1/public/yifeng/Github/anal/'+cancer+'_baseline1.txt', 'w')
    k = 0
    for sga1 in sga2deg.keys():
        #traversedset.add(sga1)
        for sga2 in sga2deg.keys():
            #if sga2 != sga1:
            #    if sga2 in traversedset: continue 
            if sga1 == sga2: continue

            overlap = len(sga2deg[sga1].intersection(sga2deg[sga2]))
            if overlap == 0: continue
            #if overlap == 0: continue
            k += 1

            print >> f, sga1+'\t'+sga2+'\t'+str(overlap)
            #sga2sgaAaff[sga1].add((sga2, overlap))
        #degset = sga2deg[sga]
    numsga = len(sga2deg.keys())
    print k
    #print numsga, numsga*numsga, k
    #print len(sga2sgaAaff), len(sga2sgaAaff.keys())

    f.close()
    #return sga2sgaAaff

if __name__ == '__main__':

    for cancer in ['pancan', 'brca', 'gbm', 'ov']:
    #for cancer in ['brca']:
        getAffinityBaseline1(cancer)

#EOF.
