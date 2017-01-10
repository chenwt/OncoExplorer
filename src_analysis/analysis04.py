from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    '''
    Calculating Jaccard weight.
    python analysis04.py --inputData /usr1/public/yifeng/Github/outputData --outputData /usr1/public/yifeng/Github/outputData --cancer brca
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
    tmp_graph = dd()

    print 'reading from... {}'.format(path_inputData+'/ensemble.txt')
    k = 0
    with open(path_inputData+'/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
            if can != cancer: continue
            k = k+1
            #print patid, can, sga, deg, prob
            if k%100000 == 0: print k 
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            else:
                tmp_graph[(sga,deg)] += 1.0
    print '#readin_edge = {}'.format(k)

    f = open(path_outputData+'/graph_'+cancer+'.txt', 'w')
    print >> f, 'SGA\tDEG\tWeight'
    for line in tmp_graph.keys():
        sga,deg,weight = line[0],line[1],tmp_graph[(line[0],line[1])]
        print >> f, sga+'\t'+deg+'\t'+str(weight)
    f.close()



    sga2deg = dd(set)

    sgaWsga = dd(float)

    for line in tmp_graph.keys():
        sga,deg = line[0],line[1]
        sga2deg[sga].add(deg)

    # count reciprocal distance matrix
    checked = set()
    for sga1 in sga2deg.keys():
        for sga2 in sga2deg.keys():
            if sga2 == sga1: continue
            if (sga1, sga2) in checked: continue
            checked.add((sga1,sga2))
            checked.add((sga2,sga1))

            if len(sga2deg[sga1].intersection(sga2deg[sga2])) == 0: continue
            # Jaccard weight/coefficient
            sgaWsga[(sga1,sga2)] = 1.0* len(sga2deg[sga1].intersection(sga2deg[sga2])) / len(sga2deg[sga1].union(sga2deg[sga2]))


    print 'len(link)=%d'%len(sgaWsga)

    gene_set = set()
    for line in sgaWsga.keys():
        gene_set.add(line[0])
        gene_set.add(line[1])

    print 'len(node)=%d'%len(gene_set)

    f = open(path_outputData+'/linked_'+cancer+'.txt', 'w')
    print >> f, 'Source\tTarget\tType\tWeight'
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        wt = sgaWsga[(sga1,sga2)]
        print >> f, sga1+'\t'+sga2+'\tUndirected\t'+str(sgaWsga[(sga1,sga2)])
    f.close()    


    gene2id = dd()
    k = 0
    for line in gene_set:
        k += 1
        gene2id[line] = k
    print 'k = %d'%k

    f = open(path_outputData+'/linked_'+cancer+'_n.txt', 'w')
    print >> f, 'Source\tTarget\tWeight'
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        wt = sgaWsga[(sga1,sga2)]
        print >> f, str(gene2id[sga1])+'\t'+str(gene2id[sga2])+'\t'+str(sgaWsga[(sga1,sga2)])
    f.close()    

    f = open(path_outputData+'/lut_'+cancer+'.txt', 'w')
    print >> f, 'Gene\tId'
    for line in gene_set:
        print >> f, line+'\t'+str(gene2id[line])
    f.close()    
#EOF.
