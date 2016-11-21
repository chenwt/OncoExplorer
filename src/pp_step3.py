from collections import defaultdict as dd
import io
import os
import random
import argparse

def writeSample_Proppr(path, filename, relname, src2dst_list, dst_corpus):
    # print 'saving to {}...'.format(path+'/'+filename)

    f = io.open(path+'/tmp','w')
    for itr, src in enumerate(src2dst_list):
        if itr%1000 == 0:
            print itr
        dst = src2dst_list[src]
        f.write( u'%s(%s,Y)'%(relname,src) )
        for gene in dst_corpus:
            if gene in dst:
                f.write(u'\t+')
            else:
                f.write(u'\t-')
            f.write( u'%s(%s,%s)'%(relname,src,dst) )
        f.write(u'\n')
    f.close()

    examples = []
    for line in open(path+'/tmp', 'r'):
        line = line.strip()
        examples.append(line)

    print 'len(examples) = {}'.format(len(examples))

    SEED = 6
    random.seed(SEED)
    random.shuffle(examples)
    os.remove(path+'/tmp');

    f = open(path+'/'+filename, 'w')
    for line in examples:
        print >> f, line
    f.close()

def writeSample_Tensorlog(path, filename, relname, src2dst_list):
    f = io.open(path+'/'+filename,'w')
    for src in src2dst_list.keys():
        dst = src2dst_list[src]
        f.write(u'%s\t%s'%(relname,src))
        for gene in dst:
            f.write(u'\t%s'%gene)
        f.write(u'\n')

if __name__ == '__main__':
    ''' Third step to generate the files for ProPPR and TensorLog.
    Usage: python pp_step3.py --inputData <dir contatins train/test/graph.txt> --TensorlogData <dir to contain input of TensorLog> --PropprData <dir to contain input of ProPPR> 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData',     help = 'directory that contatins train/test/graph/SGA/DEG.txt', type = str)
    parser.add_argument('--TensorlogData', help = 'directory to contain the input data of TensorLog', type = str)
    parser.add_argument('--PropprData',    help = 'directory to contain the input data of ProPPR', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_TensorlogData = args.TensorlogData
    path_PropprData = args.PropprData

    sga_corpus = set()
    with open(path_inputData+'/SGA.txt') as f:
        for line in f:
            line = line.strip()
            sga_corpus.add(line)

    deg_corpus = set()
    with open(path_inputData+'/DEG.txt') as f:
        for line in f:
            line = line.strip()
            deg_corpus.add(line)

    gene_corpus = set()
    with open(path_inputData+'/Gene.txt') as f:
        for line in f:
            line = line.strip()
            gene_corpus.add(line)

    f = open(path_PropprData+'/isSGA.cfacts', 'w')
    for line in sga_corpus:
        print >> f, 'isSGA\t'+line
    f.close()

    f = open(path_PropprData+'/isDEG.cfacts', 'w')
    for line in deg_corpus:
        print >> f, 'isDEG\t'+line
    f.close()

    train_sga2deg_list = dd(list)
    train_deg2sga_list = dd(list)

    f = open(path_inputData+'/train.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        sga,deg = line[0], line[1]
        train_sga2deg_list[sga].append(deg)
        train_deg2sga_list[deg].append(sga)
    f.close()

    test_sga2deg_list = dd(list)
    test_deg2sga_list = dd(list)

    f = open(path_inputData+'/test.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        sga,deg = line[0], line[1]
        test_sga2deg_list[sga].append(deg)
        test_deg2sga_list[deg].append(sga)
    f.close()

    #writeSample_Proppr(path_PropprData, 'train_causes.examples', 'causes', train_sga2deg_list, deg_corpus)
    #writeSample_Proppr(path_PropprData, 'test_causes.examples', 'causes', test_sga2deg_list, deg_corpus)

    #writeSample_Proppr(path_PropprData, 'train_causedBy.examples', 'causedBy', train_deg2sga_list, sga_corpus)
    #writeSample_Proppr(path_PropprData, 'test_causedBy.examples', 'causedBy', test_deg2sga_list, sga_corpus)


    
    writeSample_Tensorlog(path_TensorlogData,'train_causedBy.exam','causedBy',train_deg2sga_list)
    writeSample_Tensorlog(path_TensorlogData,'test_causedBy.exam','causedBy',test_deg2sga_list)


    graph = list()
    f = open(path_inputData+'/graph.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        sga,deg = line[0], line[1]
        graph.append((sga,deg))
    f.close()    

    f = open(path_TensorlogData+'/pathway_causedBy.cfacts', 'w')
    for line in graph:
        print >> f, 'drivenBy\t'+line[1]+'\t'+line[0]
    print >> f, ''
    for line in sga_corpus:
        print >> f, 'isSAG\t'+line
    print >> f, ''
    for deg in train_deg2sga_list.keys():
        for sga in train_deg2sga_list[deg]:
            print >> f, 'train\t'+deg+'\t'+sga
    print >> f, ''
    for deg in test_deg2sga_list.keys():
        for sga in test_deg2sga_list[deg]:
            print >> f, 'test\t'+deg+'\t'+sga
    print >> f, ''
    print >> f, 'rule\tf1'
    print >> f, 'rule\tf2'
    print >> f, ''
    for gene in gene_corpus:
        print >> f, 'src\t'+gene
    print >> f, ''
    for gene in gene_corpus:
        print >> f, 'dst\t'+gene
    f.close()
    print 'Done!'








#     # (sga,deg) -> prob
#     #sga2deg2prob_train = dd()
#     # set of existing (sga,deg)
#     set_train, set_test, set_graph = set(), set(), set()
#     # \Pi (1-prob_i) of specific (sga,deg)
#     tmp_train, tmp_test, tmp_graph = dd(), dd(), dd()
#     # \Sum prob_i of specific (sga,deg)
#     #sum_train, sum_test, sum_graph = dd(), dd(), dd()    

#     print 'reading from... {}'.format(path_inputData+'/ensemble.txt')
#     k = 0
#     with open(path_inputData+'/ensemble.txt') as f:
#         next(f)
#         for line in f:
#             line = line.strip().split('\t')
#             patid, sga, deg, prob = int(line[1]),line[4],line[5],line[6]
#             if float(prob) < threshold_readin: continue
#             k = k+1
#             if k%1000000 == 0: print k 
#             if patid in patid_train:
#                 if (sga,deg) not in set_train:
#                     tmp_train[(sga,deg)] = 1.0
#                     set_train.add((sga,deg))
#                 tmp_train[(sga,deg)] *= (1.0 - float(prob))
#                 #sum_train[(sga,deg)] += float(prob)
#             elif patid in patid_test:
#                 if (sga,deg) not in set_test:
#                     tmp_test[(sga,deg)] = 1.0
#                     set_test.add((sga,deg))
#                 tmp_test[(sga,deg)] *= (1.0 - float(prob))
#                 #sum_test[(sga,deg)] += float(prob)
#             else: # patid in patid_graph
#                 if (sga,deg) not in set_graph:
#                     tmp_graph[(sga,deg)] = 1.0
#                     set_graph.add((sga,deg))
#                 tmp_graph[(sga,deg)] *= (1.0 - float(prob))
#                 #sum_graph[(sga,deg)] += float(prob)
#     print '#readin_edge = {}'.format(k)

#     f = open(path_outputData+'/train.txt', 'w')
#     for line in tmp_train.keys():
# #        print line
#         sga,deg,prob = line[0],line[1],1.0-tmp_train[(line[0],line[1])]
#         #print sga,deg,str(prob)
#         if prob > threshold_noisyor:
#             print >> f, sga+'\t'+deg+'\t'+str(prob)
#     f.close()

#     f = open(path_outputData+'/test.txt', 'w')
#     for line in tmp_test.keys():
#         sga,deg,prob = line[0],line[1],1.0-tmp_test[(line[0],line[1])]
#         if prob > threshold_noisyor:
#             print >> f, sga+'\t'+deg+'\t'+str(prob)
#     f.close()

#     f = open(path_outputData+'/graph.txt', 'w')
#     for line in tmp_graph.keys():
#         sga,deg,prob = line[0],line[1],1.0-tmp_graph[(line[0],line[1])]
#         if prob > threshold_noisyor:
#             print >> f, sga+'\t'+deg+'\t'+str(prob)
#     f.close()
#End
