from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    '''
    python analysis02.py --inputData /usr1/public/yifeng/Github/outputData --outputData /usr1/public/yifeng/Github/outputData --cancer brca
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    parser.add_argument('--cancer', help = 'abbr of cancer name', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    cancer = args.cancer

    # set of existing (sga,deg)
    set_graph = set()
    # \Pi (1-prob_i) of specific (sga,deg)
    tmp_graph = dd()

    print 'reading from... {}'.format(path_inputData+'/ensemble.txt')
    k = 0
    with open(path_inputData+'/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
            #if float(prob) < threshold_readin: continue
            if can != cancer: continue
            k = k+1
            #print patid, can, sga, deg, prob
            if k%100000 == 0: print k 
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            else:
                tmp_graph[(sga,deg)] += 1.0
            #tmp_graph[(sga,deg)] *= (1.0 - float(prob))
                #sum_graph[(sga,deg)] += float(prob)
    print '#readin_edge = {}'.format(k)

    f = open(path_outputData+'/graph_'+cancer+'.txt', 'w')
    print >> f, 'SGA\tDEG\tWeight'
    for line in tmp_graph.keys():
        sga,deg,weight = line[0],line[1],tmp_graph[(line[0],line[1])]
        #sga,deg,prob = line[0],line[1],1.0-tmp_graph[(line[0],line[1])]
        #if prob > threshold_noisyor:
        print >> f, sga+'\t'+deg+'\t'+str(weight)
    f.close()



    sga2deg = dd(list)

    deg2sga = dd(list)

    sgaWsga = dd(int)

    for line in tmp_graph.keys():
        sga,deg = line[0],line[1]
        sga2deg[sga].append(deg)
        deg2sga[deg].append(sga)

    # count reciprocal distance matrix
    for sga in sga2deg.keys():
        for deg in sga2deg[sga]:
            for sga2 in deg2sga[deg]:
                sgaWsga[(sga,sga2)] += 1.0
    # remove digonal elements
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        if sga1 == sga2: del sgaWsga[(sga1,sga2)]
    #print len(sgaWsga)    

    # remove half of matrix
    del_set = set()
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        if (sga1, sga2) not in del_set:
            del_set.add((sga2, sga1))
            del sgaWsga[(sga1,sga2)]
    print 'len(link)=%d'%len(sgaWsga)

    # for line in sgaWsga.keys():
    #     if sgaWsga[line] < threshold: del sgaWsga[line]
    # print 'len(flink)=%d'%len(sgaWsga)

    gene_set = set()
    for line in sgaWsga.keys():
        gene_set.add(line[0])
        gene_set.add(line[1])

    print 'len(node)=%d'%len(gene_set)

    f = open(path_outputData+'/linked_'+cancer+'.txt', 'w')
    print >> f, 'Source\tTarget\tType\tWeight'
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        scale = 100.0
        print >> f, sga1+'\t'+sga2+'\tUndirected\t'+str(sgaWsga[(sga1,sga2)]/scale)
    f.close()    

#End
