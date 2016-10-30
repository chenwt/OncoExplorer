from collections import defaultdict as dd
import io
import os
import random

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

# def get_sga2deg(path_sga2deg_all, path_sga2deg_train, patid_train):
#     print 'reading from: {}...'.format(path_sga2deg_all)
#     sga2deg = dd(float)
#     sga2deg_str = []
#     sga2deg_all_list = set()
#     sga_corpus = set()
#     deg_corpus = set()
#     for line in open(path_sga2deg_all, 'r'):
#         values = line.strip().split('\t')
#         patid,sga,deg,prob = int(values[0]),values[1],values[2],float(values[3])
#         sga2deg_all_list.add((sga,deg))
#         if patid not in patid_train: continue
#         sga2deg[(sga,deg)]+= prob
#         sga_corpus.add(sga)
#         deg_corpus.add(deg)

#     for row in sga2deg.keys():
#         sga2deg_str.append('sga2deg'+'\t'+row[0]+'\t'+row[1]+'\t'+str(sga2deg[row]))

#     print 'saving to {}...'.format(path_sga2deg_train)
#     f = open(path_sga2deg_train,'w')
#     for row in sga2deg_str:
#         print >> f, row
#     f.close

#     sga2deg = sga2deg.keys()
#     return sga2deg,sga_corpus,deg_corpus,sga2deg_all_list

def writeSample(path, filename, sga2deglist, deg_corpus):
    with io.open(path+'/tmp','w') as file:
        for itr, sga in enumerate(sga2deglist):
            if itr%1000 == 0:
                print itr
            deg = sga2deglist[sga]
            file.write(u'pathTo(%s,Y)'%sga)
            # TODO:
            for gene in deg_corpus:
                # TODO:
                if gene in deg:
                    file.write(u'\t+')
                    #i += 1
                else:
                    file.write(u'\t-')
                    #j += 1
                file.write(u'pathTo(%s,%s)'%(sga,gene))
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


def step1(path_inputData,path_outputData):
    '''
    first step of preparing data.
    '''
    print 'extracting datasets in patient-wise manner...'

    dest = path_outputData+'/pathway'
    if not os.path.exists(path_outputData):
        os.makedirs(path_outputData)
    if not os.path.exists(dest):
        os.makedirs(dest)
    if not os.path.exists(dest+'/pathway_processed'):
        os.makedirs(dest+'/pathway_processed')

    # Read raw data from dumped .sql files.
    path_tdi = path_inputData+'/TDI_Results.sql'
    path_sga = path_inputData+'/SGAs.sql'
    path_deg = path_inputData+'/DEGs.sql'
    path_gen = path_inputData+'/Genes.sql'

    patid2sgaid2degid2prob = readTDI_tuple(path_tdi,1,2,4,5)
    sgaid2genid = readSQL(path_sga,0,2)
    degid2genid = readSQL(path_deg,0,2)
    genid2gen = readSQL(path_gen,0,1)

    # isDEG.cfacts
    deg_corpus = set()
    for _, genid in degid2genid.iteritems():
        gene = genid2gen[genid][1:-1].lower()
        deg_corpus.add(gene)

    path_deg = dest+'/pathway_processed/isDEG_all.cfacts'
    print 'saving to {}...'.format(path_deg)
    f = open(path_deg,'w')
    for gene in deg_corpus:
        print >> f, 'isDEG\t'+gene
    f.close
    print 'len(deg_corpus) = {}'.format(len(deg_corpus))

    # isSGA.cfacts
    sga_corpus = set()
    for _, genid in sgaid2genid.iteritems():
        if genid == 'NULL': continue
        gene = genid2gen[genid][1:-1].lower()
        sga_corpus.add(gene)

    path_sga = dest+'/pathway_processed/isSGA_all.cfacts'
    print 'saving to {}...'.format(path_sga)
    f = open(path_sga,'w')
    for gene in sga_corpus:
        print >> f, 'isSGA\t'+gene
    f.close
    print 'len(sga_corpus) = {}'.format(len(sga_corpus))

    print 'mapping from ids to genes...'
    sga2deg_all = list()

    for row in patid2sgaid2degid2prob:
        patid = row[0]
        sgaid_tmp = sgaid2genid[row[1]]
        degid_tmp = degid2genid[row[2]]
        # prob = float(row[3])
        prob = row[3]
        # sga is unit, no corresponding genid.
        if sgaid_tmp == 'NULL': continue
        sga = genid2gen[sgaid_tmp]
        deg = genid2gen[degid_tmp]
        # normalize to lower case gene representation.
        sga = sga[1:-1].lower()
        deg = deg[1:-1].lower()

        # sga and deg should be different.
        if sga != deg:
            sga2deg_all.append((patid,sga,deg,prob))

    path_sga2deg = dest+'/pathway_processed/sga2deg_all.graph'
    save2txt(path_sga2deg, sga2deg_all)

    print 'len(sga2deg_all) = {}'.format(len(sga2deg_all))

    print 'Done step 1!'


def step2(path_outputData):

    print 'preparing training and test set...'

    root = path_outputData
    dest = root+'/pathway'


    # get the training(?) edgesself.

    NUM_PAT = 4468
    patid_list = range(1,NUM_PAT+1)
    SEED = 888
    random.seed(SEED)
    random.shuffle(patid_list)

    thresh1 = 0.1
    thresh2 = 0.2
    cut1 = int(thresh1*NUM_PAT)
    cut2 = int(thresh2*NUM_PAT)
    patid_train = patid_list[0:cut1]
    patid_test = patid_list[cut1:cut2]
    patid_remain = patid_list[cut2:]

    # check
    print len(patid_train),len(patid_test),len(patid_remain),len(patid_train)+len(patid_test)+len(patid_remain)

    path_sga2deg_all = dest+'/pathway_processed/sga2deg_all.graph'
    path_sga2deg_train = dest+'/train'
    path_sga2deg_test = dest+'/test'
    path_sga2deg_remain = dest+'/remain'


    print 'reading from: {}...'.format(path_sga2deg_all)

    sga2deg_train = dd(float)
    sga2deg_test = dd(float)
    sga2deg_remain = dd(float)

    sga2deg_remain_imp = dd(float)

    sga_train_corpus = set()
    deg_train_corpus = set()
    sga_test_corpus = set()
    deg_test_corpus = set()
    sga_remain_corpus = set()
    deg_remain_corpus = set()

    count = 0

    sga2pat = dd(set)

    for line in open(path_sga2deg_all, 'r'):
        count += 1
        if count%1000000 == 0:
            print count
        values = line.strip().split('\t')
        patid,sga,deg,prob = int(values[0]),values[1],values[2],float(values[3])
        if patid in patid_train:
            sga_train_corpus.add(sga)
            deg_train_corpus.add(deg)
            sga2deg_train[(sga,deg)] += prob
        elif patid in patid_test:
            sga_test_corpus.add(sga)
            deg_test_corpus.add(deg)
            sga2deg_test[(sga,deg)] += prob
        elif patid in patid_remain:
            sga_remain_corpus.add(sga)
            deg_remain_corpus.add(deg)
            sga2deg_remain[(sga,deg)] += prob

            sga2pat[sga].add(patid)
        else:
            print 'error!!!'


    print 'saving to {}...'.format(path_sga2deg_train)
    f = open(path_sga2deg_train, 'w')
    for row in sga2deg_train.keys():
        print >> f, row[0]+'\t'+row[1]+'\t'+str(sga2deg_train[row])
    f.close()

    print 'saving to {}...'.format(path_sga2deg_test)
    f = open(path_sga2deg_test, 'w')
    for row in sga2deg_test.keys():
        print >> f, row[0]+'\t'+row[1]+'\t'+str(sga2deg_test[row])
    f.close()

    print 'saving to {}...'.format(path_sga2deg_remain)
    f = open(path_sga2deg_remain, 'w')
    for row in sga2deg_remain.keys():
        print >> f, row[0]+'\t'+row[1]+'\t'+str(sga2deg_remain[row])
    f.close()


    print 'sga2deg_train\tsga2deg_test\tsga2deg_remain'
    print '{}\t\t{}\t\t{}'.format(len(sga2deg_train.keys()),len(sga2deg_test.keys()),len(sga2deg_remain.keys()))
    print '\t\t{}\t\t{}'.format(len(set(sga2deg_train.keys()).intersection(set(sga2deg_test.keys()))), \
        len(set(sga2deg_train.keys()).intersection(set(sga2deg_remain.keys()))))

    print 'sga_train\tsga_test\tsga_remain'
    print '{}\t\t{}\t\t{}'.format(len(sga_train_corpus),len(sga_test_corpus),len(sga_remain_corpus))
    print '\t\t{}\t\t{}'.format(len(sga_train_corpus.intersection(sga_test_corpus)), \
        len(sga_train_corpus.intersection(sga_remain_corpus)))

    print 'deg_train\tdeg_test\tdeg_remain'
    print '{}\t\t{}\t\t{}'.format(len(deg_train_corpus),len(deg_test_corpus),len(deg_remain_corpus))
    print '\t\t{}\t\t{}'.format(len(deg_train_corpus.intersection(deg_test_corpus)), \
        len(deg_train_corpus.intersection(deg_remain_corpus)))

    path_graph = dest+'/sga2deg.graph'
    print 'saving to {}...'.format(path_graph)
    f = open(path_graph,'w')
    for row in sga2deg_remain.keys():
        sga,deg,prob = row[0],row[1],str(sga2deg_remain[row])
        print >> f, 'sga2deg\t'+sga+'\t'+deg+'\t'+prob

    f.close

    # get the labels.cfacts
    path_label = dest+'/isDEG.cfacts'
    print 'saving to {}...'.format(path_label)
    f = open(path_label,'w')
    for deg in deg_train_corpus:
        print >> f, 'isDEG\t'+deg
    f.close
    print 'Done step 2!'

def step3(path_outputData):

    print 'preparing training and test set...'

    root = path_outputData
    dest = root+'/pathway'

    deg_corpus_train = set()
    path_deg_corpus_train = dest+'/isDEG.cfacts'
    print 'reading from: {}...'.format(path_deg_corpus_train)
    for line in open(path_deg_corpus_train, 'r'):
        values = line.strip().split('\t')
        deg_corpus_train.add(values[1])

    sga2deglist_train = dd(list)
    path_train = dest+'/train'
    print 'reading from: {}...'.format(path_train)
    for line in open(path_train, 'r'):
        values = line.strip().split('\t')
        sga,deg,prob = values[0],values[1],float(values[2])
        sga2deglist_train[sga].append(deg)

    writeSample(dest, 'train.examples', sga2deglist_train, deg_corpus_train)

    # generate the testing set.
    deg_corpus_test = set()
    path_deg_corpus_test = dest+'/isDEG.cfacts'

    sga2deglist_test = dd(list)
    path_test= dest+'/test'
    print 'reading from: {}...'.format(path_test)
    for line in open(path_test, 'r'):
        values = line.strip().split('\t')
        sga,deg,prob = values[0],values[1],float(values[2])
        sga2deglist_test[sga].append(deg)
        deg_corpus_test.add(deg)

    writeSample(dest, 'test.examples', sga2deglist_test, deg_corpus_test)

    print 'Done step 3!'

if __name__ == '__main__':

    root = '/home/yifengt/Github'
    path_inputData = root+'/inputData'
    path_outputData = root+'/outputData'

    # step1(path_inputData,path_outputData)
    # step2(path_outputData)
    step3(path_outputData)





    



