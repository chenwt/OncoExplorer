from collections import defaultdict as dd
import io
import os
import random
import argparse

def readTDI_tuple(path, pos_patient, pos_sga, pos_deg, pos_prob):
    '''
    Collect all the (patient_id, sgaid, degid) tuples.
    path: dir of .sql file.
    pos_patient: position of patid.
    pos_sga: position of sgaid.
    pos_deg: position of degid.
    '''
    print 'reading from: {}...'.format(path)
    sga2deg = []
    for line in open(path, 'r'):
        values = line.split('),(')
        if 'INSERT INTO' in values[0]:
            values[-1] = values[-1][:-3]
            tmp = values[0].split('(')
            values[0] = tmp[1]
            for val in values:
                row = val.split(',')
                if row[pos_sga] != 'NULL':
                    sga2deg.append((row[pos_patient],row[pos_sga],row[pos_deg],row[pos_prob]))
    print 'len(patid2sgaid2degid2prob) = {}'.format(len(sga2deg))
    return sga2deg

def readSQL(path, pos_src, pos_dist):
    '''
    Read LUT from .sql file.
    src2dist: a dictionary
    '''
    print 'reading from: {}...'.format(path)
    src2dist = {}
    for line in open(path, 'r'):
        values = line.split('),(')
        if 'INSERT INTO' in values[0]:
            # modify first and last string in list.
            values[-1] = values[-1][:-3]
            tmp = values[0].split('(')
            values[0] = tmp[1]
            for val in values:
                row = val.split(',')
                src2dist[row[pos_src]] = row[pos_dist]

    print 'len = {}'.format(len(src2dist))
    return src2dist

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

def save2txt_list(path, table):
    '''
    path: the filename to be saved.
    table: list/set of string, data to be saved.
    '''
    print 'saving to {}...'.format(path)
    f = open(path, 'w')
    for row in table:
        print >> f, row
    f.close()

def writeSample(path, filename, sga2deglist, deg_corpus):
    with io.open(path+'/tmp','w') as file:
        for itr, sga in enumerate(sga2deglist):
            if itr%1000 == 0:
                print itr
            deg = sga2deglist[sga]
            file.write(u'causedBy(%s,Y)'%sga)
            # TODO:
            for gene in deg_corpus:
                # TODO:
                if gene in deg:
                    file.write(u'\t+')
                    #i += 1
                else:
                    file.write(u'\t-')
                    #j += 1
                file.write(u'causedBy(%s,%s)'%(sga,gene))
            file.write(u'\n')

    examples = []
    for line in open(path+'/tmp', 'r'):
        line = line.strip()
        examples.append(line)

    print 'len(samples) = {}'.format(len(examples))

    SEED = 666
    random.seed(SEED)
    random.shuffle(examples)
    os.remove(path+'/tmp');
    path_out = path+'/'+filename
    save2txt_list(path_out,examples)   


if __name__ == '__main__':
    

    # root = '/home/yifengt/Github'
    # path_inputData = root+'/inputData'
    # path_outputData = root+'/outputData'

    path_inputData = '/remote/curtis/yifengt/outputData/pathway_wwc'
    path_outputData = path_inputData


    sga2deg = dd(str)
    for line in open(path_inputData+'/drives.cfacts', 'r'):
        line = line.strip().split('\t')
        #print line
        sga,deg,prob = line[1],line[2],line[3]
        sga2deg[(sga,deg)] = prob
        #print sga,deg,prob

    f = open(path_outputData+'/deg2sga.graph', 'w')
    for row in sga2deg.keys():
        sga,deg = row[0],row[1]
        prob = sga2deg[(sga,deg)]
        print >> f, 'deg2sga'+'\t'+deg.lower()+'\t'+sga.lower()+'\t'+prob
#        print prob
    f.close()


    sga_corpus_test = set()
    deg2sgalist_test = dd(list)
    path_test = path_inputData+'/tdi-test.exam'
    print 'reading from: {}...'.format(path_test)
    for line in open(path_test, 'r'):
        values = line.strip().split('\t')
        deg = values[1]
        for i in range(len(values)):
            if i <= 1: continue
            sga = values[i]
            sga_corpus_test.add(sga.lower())
            deg2sgalist_test[deg.lower()].append(sga.lower())


    print len(sga_corpus_test)
    path_isSGA = path_outputData+'/isSGA.cfacts'
    f = open(path_isSGA,'w')
    for sga in sga_corpus_test:
        print >> f, 'isSGA\t'+sga
    f.close

    path_test0 = path_outputData+'/test'
    f = open(path_test0, 'w')
    for deg in deg2sgalist_test.keys():
        sgalist = deg2sgalist_test[deg]
        for sga in sgalist:
            print >> f, deg+'\t'+sga
    f.close

    writeSample(path_outputData, 'test.examples', deg2sgalist_test, sga_corpus_test)

# Usage:
# python prepare_pathway_clean.py --inputData '/remote/curtis/yifengt/inputData' --outputData '/remote/curtis/yifengt/outputData'

