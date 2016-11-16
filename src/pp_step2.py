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
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    noisyor = args.noisyor
    threshold = args.threshold


    # shuffle training and test

    print 'reading from... {}'.format(path_inputData+'/ensemble.txt')
    for line in open(path_inputData+'/ensemble.txt', 'r'):
        line = line.strip().split('\t')
        pat, patid, can, canid, sga, deg, prob = line


#End
