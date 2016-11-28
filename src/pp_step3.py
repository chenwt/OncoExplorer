from collections import defaultdict as dd
import io
import os
import random
import argparse

def writetestSample_Proppr(path, filename, relname, src_corpus):
    f = io.open(path+'/'+filename,'w')
    for src in sga_corpus:
        f.write( u'%s(%s,Y)\n'%(relname,src) )
    f.close()

def writeSample_Tensorlog(path, filename, relname, src2dst_list):
    f = io.open(path+'/'+filename,'w')
    for src in src2dst_list.keys():
        dst = src2dst_list[src]
        f.write(u'%s\t%s'%(relname,src))
        for gene in dst:
            f.write(u'\t%s'%gene)
        f.write(u'\n')

if __name__ == '__main__':
    ''' Third step to generate the files for ProPPR and TensorLog.
    Usage: python pp_step3.py --inputData <dir contatins train/test/graph.txt> --TensorlogData <dir to contain input of TensorLog> --PropprData <dir to contain input of ProPPR> 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData',     help = 'directory that contatins train/test/graph/SGA/DEG.txt', type = str)
    parser.add_argument('--TensorlogData', help = 'directory to contain the input data of TensorLog', type = str)
    parser.add_argument('--PropprData',    help = 'directory to contain the input data of ProPPR', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_TensorlogData = args.TensorlogData
    path_PropprData = args.PropprData

    sga_corpus = set()
    with open(path_inputData+'/SGA.txt') as f:
        for line in f:
            line = line.strip()
            sga_corpus.add(line)

    writetestSample_Proppr(path_PropprData, 'test_linked.examples', 'linked', sga_corpus)

    graph = list()
    f = open(path_inputData+'/graph.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        sga,deg = line[0], line[1]
        graph.append((sga,deg))
    f.close()    

    f = open(path_PropprData+'/pathway.graph', 'w')
    for line in graph:
        sga,deg = line[0], line[1]
        print >> f, 'drives\t'+sga+'\t'+deg
    f.close()
    f = open(path_PropprData+'/pathwayr.graph', 'w')
    for line in graph:
        sga,deg = line[0], line[1]
        print >> f, 'drivenBy\t'+deg+'\t'+sga
    f.close()

    print 'Done!'
#EOF.
