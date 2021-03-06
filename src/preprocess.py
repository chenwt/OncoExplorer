from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':

    # root = '/home/yifengt/Github'
    # path_inputData = root+'/inputData'
    # path_outputData = root+'/outputData'
    root = '/home/yifengt/Github'
    path_inputData = root+'/inputData'
    path_outputData = root+'/outputData'

    pat2sga2deg = set()
#TODO: first line to be deleted.
    print 'reading into filtered results...'
    for line in open(path_inputData+'/TDI_Results_filter.csv', 'r'):
        line = line.strip().split(',')
        pat,sga,deg = line[2],line[0],line[1]
        pat2sga2deg.add((pat,sga,deg))
#        print pat, sga, deg
    pat2sga2deg2prob = dd(str)
    for line in open(path_inputData+'/TDI_Results_new.csv', 'r'):
        line = line.strip().split(',')
        pat,sga,deg,prob = line[0],line[1],line[2],line[3]
        if (pat,sga,deg) in pat2sga2deg:
            pat2sga2deg2prob[(pat,sga,deg)] = prob
#        print pat,sga,deg,prob

    f = open(path_outputData+'/pat2deg2sga2prob.txt', 'w')
    for line in pat2sga2deg2prob.keys():
        pat,sga,deg,prob = line[0],line[1],line[2],pat2sga2deg2prob[line]
        print >> f, pat+'\t'+sga+'\t'+deg+'\t'+prob 
    f.close()


# Usage:
# python prepare_pathway_clean.py --inputData '/remote/curtis/yifengt/inputData' --outputData '/remote/curtis/yifengt/outputData'

