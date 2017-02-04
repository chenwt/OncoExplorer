# Count the degrees of SGAs and output them into the files.
from collections import defaultdict as dd
import numpy as np
import os
from scipy.io import savemat
from scipy.sparse import csr_matrix

if __name__ == '__main__':
    sga2deg = dd(set)
    #tmp_graph = dd(float)
    #set_deg = set()
    #deg2sga = dd(set)
    sga2sga = dd(set)
    k = 0
    flag = True
    with open('/usr1/public/yifeng/Github/inputData/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
            k = k+1
            sga2deg[sga].add(deg)
            #else:
                #tmp_graph[(sga,deg)] += 1.0
    set_sga = list(sga2deg.keys())

    f = open('/usr1/public/yifeng/Github/OncoExplorer/src/degree_SGA_DEG.txt', 'w')
    print >> f, 'SGA\tDegree'
    for sga in set_sga:
        print >> f, sga+'\t'+str(len(sga2deg[sga]))
    f.close()


    for sga1 in set_sga:
        for sga2 in set_sga:
            if len(sga2deg[sga1].intersection(sga2deg[sga2])) > 0:
                sga2sga[sga1].add(sga2)
                sga2sga[sga2].add(sga1)
    f = open('/usr1/public/yifeng/Github/OncoExplorer/src/degree_SGA_SGA.txt', 'w')
    print >> f, 'SGA\tDegree'
    for sga in set_sga:
        print >> f, sga+'\t'+str(len(sga2sga[sga]))
        #print >> f, sga+'\t'+str(len(sga2deg[sga]))

    f.close()

#EOF.
