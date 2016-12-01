# Extract positive grounded examples.
from collections import defaultdict as dd
import io
import os
import random
import argparse

def extract(path, pos_start):
    link = dd(float)

    for line in open(path, 'r'):
        line = line.strip().split('\t')
        if '#' in line[0]: continue
        prob = line[1]
        # get the last 'linked(,)' and strip 'linked(' and ')'
        value = line[-1][pos_start:-2]
        value = value.split(',')
        src, dst = value[0], value[1]

#        if src == dst: continue
        if float(prob) >= 0.01:
            link[(src,dst)] = float(prob)

#        link[(src,dst)] += 1
#        link[(dst,src)] += 1
    return link

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


    pi3k = set()

    f = open('/home/yifengt/Github/inputData/pi3k.txt','r')
    for line in f:
        line = line.strip().lower() 
        pi3k.add(line)
    f.close()
    pi3k_sel = dd(int) 

    weight = dd(float)
    link = extract(path_inputData+'/test_linked.solutions.txt', pos_start)
    print len(link)
    for line in link.keys():
        #if link[line] == 1: continue
        src, dst, prob = line[0],line[1],link[line]
        #if (dst,src) in link_set: continue
        #link_set.add((src,dst))
        pi3k_sel[src], pi3k_sel[dst] = 0, 0
        weight[(src, dst)] = prob
        if src in pi3k: pi3k_sel[src] = 1
        if dst in pi3k: pi3k_sel[dst] = 1
    print len(weight)
    f = open(path_outputData+'/linked.txt', 'w')
    print >> f, 'Source\tTarget\tType\tWeight'
    for line in weight.keys():
        print >> f, line[0]+'\t'+line[1]+'\tUndirected\t'+str(weight[line])
    f.close()
    f = open(path_outputData+'/node.txt', 'w')
    print >> f, 'Id\tTP53pathway'
    for line in pi3k_sel.keys():
        print >> f, line + '\t'+ str(pi3k_sel[line])
    f.close()

    print 'Done!'
#EOF.
