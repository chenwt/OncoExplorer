from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    #parser.add_argument('--threshold_readin', help = 'threshold of whether to read into a raw record of (sga,deg) pair', type = float)
    #parser.add_argument('--threshold_noisyor', help = 'threshold of whether a (sga,deg) pair exists', type = float)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    #threshold_readin = args.threshold_readin
    #threshold_noisyor = args.threshold_noisyor

    # sga : (patid, deg)
    sga2deg = dd(list)
    # (patid, deg) : sga
    deg2sga = dd(list)

    sgaWdeg = dd(int)


    print 'reading from... {}'.format(path_inputData+'/ensemble.txt')

    f = open(path_inputData+'/ensemble.txt', 'r')
    next(f)
    k = 0
    for line in f:
        line = line.strip().split('\t')
        patid, sga, deg = line[1], line[4], line[5]
        sga2deg[sga].append((patid, deg))
        deg2sga[(patid, deg)].append(sga)
        # print patid, sga, deg
        # patid, sga, deg, prob = int(line[1]),line[4],line[5],line[6]
        # if float(prob) < threshold_readin: continue
        k = k+1
        if k%100000 == 0: print k 
        # if (sga,deg) not in set_graph:
        #     tmp_graph[(sga,deg)] = 1.0
        #     set_graph.add((sga,deg))


    f.close()
    print len(sga2deg), len(deg2sga)




    # k = 0
    # with open(path_inputData+'/ensemble.txt') as f:
    #     next(f)
    #     for line in f:
    #         line = line.strip().split('\t')
    #         patid, sga, deg, prob = int(line[1]),line[4],line[5],line[6]
    #         if float(prob) < threshold_readin: continue
    #         k = k+1
    #         if k%1000000 == 0: print k 
    #         if (sga,deg) not in set_graph:
    #             tmp_graph[(sga,deg)] = 1.0
    #             set_graph.add((sga,deg))
    #         tmp_graph[(sga,deg)] *= (1.0 - float(prob))
    #             #sum_graph[(sga,deg)] += float(prob)
    # print '#readin_edge = {}'.format(k)

    # f = open(path_outputData+'/graph.txt', 'w')
    # for line in tmp_graph.keys():
    #     sga,deg,prob = line[0],line[1],1.0-tmp_graph[(line[0],line[1])]
    #     if prob > threshold_noisyor:
    #         print >> f, sga+'\t'+deg+'\t'+str(prob)
    # f.close()
#EOF.
