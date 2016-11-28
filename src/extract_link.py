# Extract positive grounded examples.
from collections import defaultdict as dd
import io
import os
import random
import argparse

def extract(path, pos_start):
    #print 'reading from: {}...'.format(path)
    link = dd(int)


    for line in open(path, 'r'):
        line = line.strip().split('\t')
        if '#' in line[0]: continue
        prob = line[1]
        # get the last 'linked(,)' and strip 'linked(' and ')'
        value = line[-1][pos_start:-2]
        value = value.split(',')
        src, dst = value[0], value[1]

#        print src, dst
        if src == dst: continue
        link[(src,dst)] += 1
        link[(dst,src)] += 1
#        if link[(dst,src)] == 1:
#            link[(dst, src)] += 1
#        else:
#            link[(src,dst)] += 1

    # src_tmp,dst_tmp = 'src_tmp','dst_tmp'
    # prob_tmp = 0.0
    # flag = 0
    # for line in open(path, 'r'):
    #     values = line.strip().split('\t')
    #     #print values
    #     if '#' in values[0]:
    #         flag = 1
    #         continue
    #     prob = values[1]
    #     val = values[-1][:-1]
    #     val = val[pos_start:-1]
    #     val = val.split(',')
    #     src,dst = val[0], val[1]
    #     prob = float(prob)

    #     T = k+1
    #     if flag == 1:
    #         src_tmp,dst_tmp,prob_tmp = src,dst,prob
    #         src2dst.add((src,dst))
    #         flag = flag + 1
    #     elif flag <= T-1: # possiblity1: a new src; possibility2: continued src
    #         # obviously, src == src_tmp:
    #         src_tmp,dst_tmp,prob_tmp = src,dst,prob
    #         src2dst.add((src,dst))
    #         # else:
    #         # do nothing
    #         flag = flag + 1
    #     elif flag == T:
    #         continue
    # print src2dst
    # print 'len(src2dst) = {}'.format(len(src2dst))
    return link



def save2txt(path, table):
    '''
    path: the filename to be saved.
    table: list of string list / set of tuples, data to be saved.
    '''
    print 'saving to {}...'.format(path)
    f = open(path, 'w')
    for row in table:
        print >> f, '\t'.join(row)
    f.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of inputData', type = str)
    parser.add_argument('--outputData', help = 'directory of outputData', type = str)
    args = parser.parse_args()
    path_inputData = args.inputData
    path_outputData = args.outputData

    print 'extracting...'

    relation = 'linked'
    pos_start = len(relation) + 1

    link = extract(path_inputData+'/test_linked.solutions.txt', pos_start)

    print 'Done!'
#EOF.
