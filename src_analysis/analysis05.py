from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    '''
    Check results.
    python analysis05.py --inputData /usr1/public/yifeng/Github/inputData --outputData /usr1/public/yifeng/Github/inputData
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData

    # set of existing (sga,deg)
    set_mine = set()
    set_old = set()
    set_new = set()

    pat_mine = set()
    pat_old = set()
    pat_new = set()
    with open(path_inputData+'/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            pat, patid, can, sga, deg, prob = line[0], int(line[1]),line[2],line[4],line[5],line[6]
            set_mine.add((pat,sga,deg))
            pat_mine.add(pat)
    f.close()

    with open(path_inputData+'/YX_TDI_SGA_DEG_Tumor.csv') as f:
        next(f)
        for line in f:
            line = line.strip().split(',')
            pat, sga, deg = line[0].lower(), line[1].lower(), line[2].lower()
            set_old.add((pat,sga,deg))
            pat_old.add(pat)
    f.close()

    with open(path_inputData+'/YX_TDIresult_PanCancerAtlas_CT01_114.csv') as f:
        next(f)
        for line in f:
            line = line.strip().split(',')
            pat, sga, deg = line[0].lower(), line[1].lower(), line[2].lower()
            set_new.add((pat,sga,deg))
            pat_new.add(pat)
    f.close()

    print '#mine={}\t#old={}\t#new={}'.format(len(set_mine), len(set_old), len(set_new))
    print '# old intersection mine ={}'.format(len(set_mine.intersection(set_old)))
    print '# old intersection new ={}'.format(len(set_old.intersection(set_new)))

    print '#mine={}\t#old={}\t#new={}'.format(len(pat_mine), len(pat_old), len(pat_new))
    print '# old intersection new ={}'.format(len(pat_old.intersection(pat_new)))
#EOF.
