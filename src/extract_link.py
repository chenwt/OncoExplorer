# Extract positive grounded examples.
from collections import defaultdict as dd
import io
import os
import random
import argparse

def extract(path, pos_start, threshold):
    link = dd(float)

    for line in open(path, 'r'):
        line = line.strip().split('\t')
        if '#' in line[0]: continue
        prob = line[1]
        # get the last 'linked(,)' and strip 'linked(' and ')'
        value = line[-1][pos_start:-2]
        value = value.split(',')
        src, dst = value[0], value[1]

        if src == dst: continue
        if float(prob) >= threshold:
            link[(src,dst)] = float(prob)
    return link

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of inputData', type = str)
    parser.add_argument('--outputData', help = 'directory of outputData', type = str)
    parser.add_argument('--threshold', help = 'threshold of taking edges', type = float)
    parser.add_argument('--pathwayName', help = 'path to the pathway file', type = str)
    args = parser.parse_args()
    path_inputData = args.inputData
    path_outputData = args.outputData
    threshold = args.threshold
    path_pathwayName = args.pathwayName
    print 'extracting...'

    relation = 'linked'
    pos_start = len(relation) + 1


    # The components of TP53 pathway.
    tp53_set = set()
    f = open(path_pathwayName,'r')
    for line in f:
        line = line.strip().lower() 
        tp53_set.add(line)
    f.close()

    # Extract the link from proppr rest.
    # threshold to decide if the weight is large enough to to take the link. 
    link = extract(path_inputData+'/test_linked.solutions.txt', pos_start, threshold)
    print len(link)

    # Get ground truth of the pathway.
    tp53_sel = dd(int) 
    for line in link.keys():
        src, dst, prob = line[0],line[1],link[line]
        tp53_sel[src], tp53_sel[dst] = 0, 0
    for line in link.keys():
        src, dst = line[0],line[1]
        if src in tp53_set: tp53_sel[src] = 1
        if dst in tp53_set: tp53_sel[dst] = 1

    


    f = open(path_outputData+'/linked.txt', 'w')
    print >> f, 'Source\tTarget\tType\tWeight'
    for line in link.keys():
        print >> f, line[0]+'\t'+line[1]+'\tDirected\t'+str(link[line])
    f.close()
    f = open(path_outputData+'/node.txt', 'w')
    print >> f, 'Id\tTP53pathway'
    for line in tp53_sel.keys():
        print >> f, line + '\t'+ str(tp53_sel[line])
    f.close()

    print 'Done!'
#EOF.
