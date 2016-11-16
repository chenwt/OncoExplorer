from collections import defaultdict as dd
import io
import os
import random
import argparse

def readSQL(path, pos_src, pos_dist):
    '''
    Read LUT from .sql file.
    src2dist: a dictionary
    '''
    print 'reading into... {}'.format(path)
    src2dist = {}
    for line in open(path, 'r'):
        values = line.split('),(')
        # reach to the line of data.
        if 'INSERT INTO' in values[0]:
            # modify first and last string in list.
            values[-1] = values[-1][:-3]
            tmp = values[0].split('(')
            values[0] = tmp[1]
            for val in values:
                row = val.split(',')
                src2dist[row[pos_src]] = row[pos_dist]

    print '#src2dist = {}'.format(len(src2dist))
    return src2dist


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData

    # patient -> (patient_id, cancer_id)
    pat2patid2canid = dd()
    print 'reading into... {}'.format(path_inputData+'/Patients.sql')
    for line in open(path_inputData+'/Patients.sql', 'r'):
        values = line.split('),(')
        # reach to the line of data.
        if 'INSERT INTO' in values[0]:
            # modify first and last string in list.
            values[-1] = values[-1][:-3]
            tmp = values[0].split('(')
            values[0] = tmp[1]
            for val in values:
                row = val.split(',')
                pat = row[1][1:-1]#remove the double quotation
                pat2patid2canid[pat] = (row[0],row[7])
    print pat2patid2canid
    print '#patients = {}'.format(len(pat2patid2canid))


    # cancer_id -> cancer
    print 'reading into... {}'.format(path_inputData+'/Cancers.sql')
    src2dist = {}
    for line in open(path, 'r'):
        values = line.split('),(')
        # reach to the line of data.
        if 'INSERT INTO' in values[0]:
            # modify first and last string in list.
            values[-1] = values[-1][:-3]
            tmp = values[0].split('(')
            values[0] = tmp[1]
            for val in values:
                row = val.split(',')
                src2dist[row[pos_src]] = row[pos_dist]

    print '#src2dist = {}'.format(len(src2dist))




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



    f = open(path_outputData+'/pat2patid2can2canid2sga2deg2prob.txt', 'w')
    for line in pat2sga2deg2prob.keys():
        pat,sga,deg,prob = line[0],line[1],line[2],pat2sga2deg2prob[line]
        print >> f, pat+'\t'+sga+'\t'+deg+'\t'+prob 
    f.close()


# Usage:
# python prepare_pathway_clean.py --inputData '/remote/curtis/yifengt/inputData' --outputData '/remote/curtis/yifengt/outputData'

