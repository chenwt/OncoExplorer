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
    parser.add_argument('--noisyor', help = 'use noisy-or or simply sum up',type = bool)
    parser.add_argument('--threshold', help = 'threshold of whether a (sga,deg) pair exists', type = float)
    parser.add_argument('--division', help = 'proportion of training patients (same as test)', type = float)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    noisyor = args.noisyor
    threshold = args.threshold
    division = args.division

    # shuffle training and test
    NUM_PAT = 4468
    patid_list = range(1,NUM_PAT+1)
    random.seed(6)
    random.shuffle(patid_list)
    cut1, cut2 = int(division*NUM_PAT), int(2*divisoin*NUM_PAT)
    patid_train, patid_test = patid_list[0:cut1], patid_list[cut1:cut2]


    # (sga,deg) -> prob
    #sga2deg2prob_train = dd()
    # set of existing (sga,deg)
    set_train, set_test, set_graph = set(), set(), set()
    # \Pi (1-prob_i) of specific (sga,deg)
    tmp_train, tmp_test, tmp_graph = dd(), dd(), dd()

    print 'reading from... {}'.format(path_inputData+'/ensemble.txt')
    for line in open(path_inputData+'/ensemble.txt', 'r'):
        line = line.strip().split('\t')
        pat, patid, can, canid, sga, deg, prob = line
        if patid in patid_train:
            if (sga,deg) not in set_train:
                tmp_train[(sga,deg)] = 1.0
                set_train.add((sga,deg))
            tmp_train[(sga,deg)] *= (1.0 - float(prob))
        elif patid in patid_test:
            if (sga,deg) not in set_test:
                tmp_test[(sga,deg)] = 1.0
                set_test.add((sga,deg))
            tmp_test[(sga,deg)] *= (1.0 - float(prob))
        else # patid in patid_graph:
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            tmp_graph[(sga,deg)] *= (1.0 - float(prob))





#End
