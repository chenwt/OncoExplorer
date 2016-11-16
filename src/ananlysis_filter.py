from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    

    # root = '/home/yifengt/Github'
    # path_inputData = root+'/inputData'
    # path_outputData = root+'/outputData'

    path_inputData = '~/Github/inputData'
    path_outputData = path_inputData

    pat2sga2deg = set()

    for line in open(path_inputData+'/TDI_Results_filter.csv', 'r'):
        line = line.strip().split(',')
        pat,sga,deg = line[2],line[0],line[1]
        pat2sga2deg.add((pat,sga,deg))

    pat2sga2deg2prob = dd(str)
    for line in open(path_inputData+'/TDI_Results_new.csv', 'r'):
        line = line.strip().split(',')
        pat,sga,deg,prob = line[0],line[1],line[3],line[4]
        if (pat,sga,deg) in pat2sga2deg:
            pat2sga2deg2prob[(pat,sga,deg)] = prob


    f = open(path_outputData+'/deg2sga.graph', 'w')
    for row in sga2deg.keys():
        sga,deg = row[0],row[1]
        prob = sga2deg[(sga,deg)]
        print >> f, 'deg2sga'+'\t'+deg.lower()+'\t'+sga.lower()+'\t'+prob
    f.close()


# Usage:
# python prepare_pathway_clean.py --inputData '/remote/curtis/yifengt/inputData' --outputData '/remote/curtis/yifengt/outputData'

