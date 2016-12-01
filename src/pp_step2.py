from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    ''' Second step to generate the graph, training and test data for ProPPR.
    Usage:
    python pp_step2.py --inputData <dir contatins ensemble.txt> --outputData <dir to contain output> 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    parser.add_argument('--threshold_readin', help = 'threshold of whether to read into a raw record of (sga,deg) pair', type = float)
    parser.add_argument('--threshold_noisyor', help = 'threshold of whether a (sga,deg) pair exists', type = float)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    threshold_readin = args.threshold_readin
    threshold_noisyor = args.threshold_noisyor

    # set of existing (sga,deg)
    set_graph = set()
    # \Pi (1-prob_i) of specific (sga,deg)
    tmp_graph = dd()

    print 'reading from... {}'.format(path_inputData+'/ensemble.txt')
    k = 0
    with open(path_inputData+'/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            patid, sga, deg, prob = int(line[1]),line[4],line[5],line[6]
            if float(prob) < threshold_readin: continue
            k = k+1
            if k%1000000 == 0: print k 
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            tmp_graph[(sga,deg)] *= (1.0 - float(prob))
                #sum_graph[(sga,deg)] += float(prob)
    print '#readin_edge = {}'.format(k)

    f = open(path_outputData+'/graph.txt', 'w')
    for line in tmp_graph.keys():
        sga,deg,prob = line[0],line[1],1.0-tmp_graph[(line[0],line[1])]
        if prob > threshold_noisyor:
            print >> f, sga+'\t'+deg+'\t'+str(prob)
    f.close()
#End
