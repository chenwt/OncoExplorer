from collections import defaultdict as dd
import io
import os
import random
import argparse


def writeSample(path, filename, sga2deglist, deg_corpus):
    # print 'saving to {}...'.format(path+'/'+filename)
    #i,j = 0,0
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

    deg_corpus = set()
    with open(path_inputData+'/DEG.txt') as f:
        for line in f:
            line = line.strip()
            print line










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
