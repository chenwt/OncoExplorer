from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData

    pat2sga2deg = set()
    #first line does not overlap, deleted automatically.
    print 'reading into...'+path_inputData+'/TDI_Results_filter_no_unit.csv'
    for line in open(path_inputData+'/TDI_Results_filter_no_unit.csv', 'r'):
        line = line.strip().split(',')
        pat,sga,deg = line[2],line[0],line[1]
        pat2sga2deg.add((pat,sga,deg))
    print '#records = {}'.format(len(pat2sga2deg))

    pat2sga2deg2prob = dd(str)
    print 'reading into...'+path_inputData+'/TDI_Results_new.csv'
    for line in open(path_inputData+'/TDI_Results_new.csv', 'r'):
        line = line.strip().split(',')
        pat,sga,deg,prob = line[0],line[1],line[2],line[3]
        if (pat,sga,deg) in pat2sga2deg:
            pat2sga2deg2prob[(pat,sga,deg)] = prob

    print '#records = {}'.format(len(pat2sga2deg2prob))
    f = open(path_outputData+'/pat2sga2deg2prob.txt', 'w')
    for line in pat2sga2deg2prob.keys():
        pat,sga,deg,prob = line[0],line[1],line[2],pat2sga2deg2prob[line]
        print >> f, pat+'\t'+sga+'\t'+deg+'\t'+prob 
    f.close()


# Usage:
# python prepare_pathway_clean.py --inputData '/remote/curtis/yifengt/inputData' --outputData '/remote/curtis/yifengt/outputData'

