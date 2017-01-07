from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    '''
    python analysis02.py --inputData /usr1/public/yifeng/Github/outputData --outputData /usr1/public/yifeng/Github/outputData --cancer brca
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    parser.add_argument('--cancer', help = 'abbr of cancer name', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    cancer = args.cancer

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
            patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
            #if float(prob) < threshold_readin: continue
            if can != cancer: continue
            k = k+1
            #print patid, can, sga, deg, prob
            if k%100000 == 0: print k 
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            else:
                tmp_graph[(sga,deg)] += 1.0
            #tmp_graph[(sga,deg)] *= (1.0 - float(prob))
                #sum_graph[(sga,deg)] += float(prob)
    print '#readin_edge = {}'.format(k)

    f = open(path_outputData+'/graph_'+cancer+'.txt', 'w')
    print >> f, 'SGA\tDEG\tWeight'
    for line in tmp_graph.keys():
        sga,deg,weight = line[0],line[1],tmp_graph[(line[0],line[1])]
        #sga,deg,prob = line[0],line[1],1.0-tmp_graph[(line[0],line[1])]
        #if prob > threshold_noisyor:
        print >> f, sga+'\t'+deg+'\t'+str(weight)
    f.close()
#End
