from collections import defaultdict as dd
import io
import os
import random
import argparse

# def writeSample_Proppr(path, filename, relname, src2dst_list, dst_corpus):
#     # print 'saving to {}...'.format(path+'/'+filename)

#     f = io.open(path+'/tmp','w')
#     for itr, src in enumerate(src2dst_list):
#         if itr%1000 == 0:
#             print itr
#         dst = src2dst_list[src]
#         f.write( u'%s(%s,Y)'%(relname,src) )
#         for gene in dst_corpus:
#             if gene in dst:
#                 f.write(u'\t+')
#             else:
#                 f.write(u'\t-')
#             f.write( u'%s(%s,%s)'%(relname,src,dst) )
#         f.write(u'\n')
#     f.close()

#     examples = []
#     for line in open(path+'/tmp', 'r'):
#         line = line.strip()
#         examples.append(line)

#     print 'len(examples) = {}'.format(len(examples))

#     SEED = 6
#     random.seed(SEED)
#     random.shuffle(examples)
#     os.remove(path+'/tmp');

#     f = open(path+'/'+filename, 'w')
#     for line in examples:
#         print >> f, line
#     f.close()
def writetestSample_Proppr(path, filename, relname, src_corpus):
    # print 'saving to {}...'.format(path+'/'+filename)

    f = io.open(path+'/'+filename,'w')
    for src in sga_corpus:
        f.write( u'%s(%s,Y)\n'%(relname,src) )
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

    # deg_corpus = set()
    # with open(path_inputData+'/DEG.txt') as f:
    #     for line in f:
    #         line = line.strip()
    #         deg_corpus.add(line)

    # gene_corpus = set()
    # with open(path_inputData+'/Gene.txt') as f:
    #     for line in f:
    #         line = line.strip()
    #         gene_corpus.add(line)

    # f = open(path_PropprData+'/isSGA.cfacts', 'w')
    # for line in sga_corpus:
    #     print >> f, 'isSGA\t'+line
    # f.close()

    # f = open(path_PropprData+'/isDEG.cfacts', 'w')
    # for line in deg_corpus:
    #     print >> f, 'isDEG\t'+line
    # f.close()

    # train_sga2deg_list = dd(list)
    # train_deg2sga_list = dd(list)

    # f = open(path_inputData+'/train.txt', 'r')
    # for line in f:
    #     line = line.strip().split('\t')
    #     sga,deg = line[0], line[1]
    #     train_sga2deg_list[sga].append(deg)
    #     train_deg2sga_list[deg].append(sga)
    # f.close()

    # test_sga2deg_list = dd(list)
    # test_deg2sga_list = dd(list)

    # f = open(path_inputData+'/test.txt', 'r')
    # for line in f:
    #     line = line.strip().split('\t')
    #     sga,deg = line[0], line[1]
    #     test_sga2deg_list[sga].append(deg)
    #     test_deg2sga_list[deg].append(sga)
    # f.close()

    #writeSample_Proppr(path_PropprData, 'train_causes.examples', 'causes', train_sga2deg_list, deg_corpus)
    writetestSample_Proppr(path_PropprData, 'test_linked.examples', 'linked', sga_corpus)






    #writeSample_Proppr(path_PropprData, 'train_causedBy.examples', 'causedBy', train_deg2sga_list, sga_corpus)
    #writeSample_Proppr(path_PropprData, 'test_causedBy.examples', 'causedBy', test_deg2sga_list, sga_corpus)


    
    # writeSample_Tensorlog(path_TensorlogData,'train.exam','causedBy',train_deg2sga_list)
    # writeSample_Tensorlog(path_TensorlogData,'test.exam','causedBy',test_deg2sga_list)


    graph = list()
    f = open(path_inputData+'/graph.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        sga,deg = line[0], line[1]
        graph.append((sga,deg))
    f.close()    

    f = open(path_PropprData+'/pathway.graph', 'w')
    for line in graph:
        sga,deg = line[0], line[1]
        print >> f, 'drives\t'+sga+'\t'+deg
    f.close()


    # f = open(path_TensorlogData+'/pathway.cfacts', 'w')
    # for line in graph:
    #     print >> f, 'drivenBy\t'+line[1]+'\t'+line[0]
    # print >> f, ''
    # for line in sga_corpus:
    #     print >> f, 'isSGA\t'+line
    # print >> f, ''
    # for deg in train_deg2sga_list.keys():
    #     for sga in train_deg2sga_list[deg]:
    #         print >> f, 'train\t'+deg+'\t'+sga
    # print >> f, ''
    # for deg in test_deg2sga_list.keys():
    #     for sga in test_deg2sga_list[deg]:
    #         print >> f, 'test\t'+deg+'\t'+sga
    # print >> f, ''
    # print >> f, 'rule\tf1'
    # print >> f, 'rule\tf2'
    # print >> f, ''
    # for gene in gene_corpus:
    #     print >> f, 'src\t'+gene
    # print >> f, ''
    # for gene in gene_corpus:
    #     print >> f, 'dst\t'+gene
    # f.close()
    print 'Done!'







#EOF.