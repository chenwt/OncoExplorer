from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    ''' 
    python analysis1.py --inputData /usr1/public/yifeng/Github/inputData --outputData /usr1/public/yifeng/Github/outputData
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData


    # patient -> (patient_id, cancer_id)
    pat2patid2canid = dd()
    print 'reading from... {}'.format(path_inputData+'/Patients.sql')
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
                pat2patid2canid[pat.lower()] = (row[0],row[7])
    print '#patients = {}'.format(len(pat2patid2canid))


    # cancer_id -> cancer
    print 'reading from... {}'.format(path_inputData+'/Cancers.sql')
    canid2can = dd()
    for line in open(path_inputData+'/Cancers.sql', 'r'):
        values = line.split('),(')
        if 'INSERT INTO' in values[0]:
            values[-1] = values[-1][:-3]
            tmp = values[0].split('(')
            values[0] = tmp[1]
            for val in values:
                row = val.split(',')
                can = row[2][1:-1]
                canid2can[row[0]] = can.lower()
    print '#cancers = {}'.format(len(canid2can))

    #filter the tdi results
    # set of (patient, sga, deg)
    pat2sga2deg = set()
    #first line does not overlap, deleted automatically.
    print 'reading from...'+path_inputData+'/TDI_Results_filter_no_unit_lc.csv'
    for line in open(path_inputData+'/TDI_Results_filter_no_unit_lc.csv', 'r'):
       line = line.strip().split(',')
       pat,sga,deg = line[2],line[0],line[1]
       if sga == 'sga_name': continue
       pat2sga2deg.add((pat,sga,deg))
    print '#records = {}'.format(len(pat2sga2deg))

    # (patient, sga, deg) -> probability
    pat2sga2deg2prob = dd(str)
    print 'reading from...'+path_inputData+'/TDI_Results_new_lc.csv'
    with open(path_inputData+'/TDI_Results_new_lc.csv') as f:
        next(f)
        for line in f:
            line = line.strip().split(',')
            pat,sga,deg,prob = line[0],line[1],line[2],line[3]
            if (pat,sga,deg) in pat2sga2deg:
                pat2sga2deg2prob[(pat,sga,deg)] = prob
    print '#records = {}'.format(len(pat2sga2deg2prob))

    # pat2sga2deg2prob = dd(str)
    # print 'reading from...'+path_inputData+'/TDI_Results_new_lc.csv'
    # with open(path_inputData+'/TDI_Results_new_lc.csv') as f:
    #     next(f)
    #     for line in f:
    #         line = line.strip().split(',')
    #         pat,sga,deg,prob = line[0],line[1],line[2],line[3]
    #         #if (pat,sga,deg) in pat2sga2deg:
    #         pat2sga2deg2prob[(pat,sga,deg)] = prob
    # print '#records = {}'.format(len(pat2sga2deg2prob))

    # for line in pat2sga2deg:
    #     if line not in pat2sga2deg2prob:
    #         print line[0]+','+line[1]+','+line[2]

    # output to ensemble files.
    print 'writing into...'+path_outputData+'/ensemble.txt'
    f = open(path_outputData+'/ensemble.txt', 'w')
    print >> f, 'patient_id\tpat_id\tcancer\tcan_id\tsga\tdeg\tprobability'
    for line in pat2sga2deg2prob.keys():
        pat,sga,deg,prob = line[0],line[1],line[2],pat2sga2deg2prob[line]
        patid, canid = pat2patid2canid[pat]
        can = canid2can[canid]

        print >> f, pat+'\t'+patid+'\t'+can+'\t'+canid+'\t'+sga+'\t'+deg+'\t'+prob 
    f.close()
#EOF.
