from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    parser.add_argument('--groundData', help = 'directory of ground truth data', type = str)
    parser.add_argument('--threshold', help = 'threshold of whether output the (sga, sga) pair', type = float)
    #parser.add_argument('--threshold_noisyor', help = 'threshold of whether a (sga,deg) pair exists', type = float)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    path_groundData = args.groundData
    threshold = args.threshold
    #threshold_noisyor = args.threshold_noisyor


    sga2deg = dd(list)

    deg2sga = dd(list)

    sgaWsga = dd(int)


    print 'reading from... {}'.format(path_inputData+'/pathway.graph')

    f = open(path_inputData+'/pathway.graph', 'r')
    k = 0
    for line in f:
        line = line.strip().split('\t')
        sga, deg = line[1], line[2]
        sga2deg[sga].append(deg)
        deg2sga[deg].append(sga)
        k = k+1
        #if k%1000 == 0: print k 

    f.close()
    #print len(sga2deg), len(deg2sga)

    # count reciprocal distance matrix
    for sga in sga2deg.keys():
        for deg in sga2deg[sga]:
            for sga2 in deg2sga[deg]:
                sgaWsga[(sga,sga2)] += 1.0
    #print len(sgaWsga)

    # remove digonal elements
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        if sga1 == sga2: del sgaWsga[(sga1,sga2)]
    print len(sgaWsga)

    # remove half of matrix
    del_set = set()
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        if (sga1, sga2) not in del_set:
            del_set.add((sga2, sga1))
            del sgaWsga[(sga1,sga2)]
    print 'len(link)=%d'%len(sgaWsga)

    for line in sgaWsga.keys():
        if sgaWsga[line] < threshold: del sgaWsga[line]
    print 'len(flink)=%d'%len(sgaWsga)


    gene_set = set()
    for line in sgaWsga.keys():
        gene_set.add(line[0])
        gene_set.add(line[1])

    print 'len(node)=%d'%len(gene_set)

    pathwayName = ['coad', 'gbm', 'blca', 'prad',
    'ucec', 'brca', 'pi3k', 'cancer']
    # name_pathway -> set of genes in KEGG pathway
    pathway_set = dd(set)
    # (name_pathway, gene) -> whether gene in name_pathway
    pathway_sel = dd(int)
    for name_pathway in pathwayName:
        f = open(path_groundData+'/'+name_pathway+'.txt','r')
        for line in f:
            line = line.strip().lower()
            pathway_set[name_pathway].add(line)
        f.close()
        for line in sgaWsga.keys():
            sga1, sga2, wt = line[0], line[1], sgaWsga[line]
            pathway_sel[(name_pathway,sga1)], pathway_sel[(name_pathway,sga2)] = 0, 0
        for line in sgaWsga.keys():
            sga1, sga2 = line[0], line[1]
            if sga1 in pathway_set[name_pathway]: pathway_sel[(name_pathway,sga1)] = 1
            if sga2 in pathway_set[name_pathway]: pathway_sel[(name_pathway,sga2)] = 1
    f = open(path_outputData+'/node.txt', 'w')
    f.write( u'Id')
    for name_pathway in pathwayName:
        f.write( u'\t%s'%name_pathway)
    f.write(u'\n')
    for sga in gene_set:
        f.write(u'%s'%(sga))
        for name_pathway in pathwayName:
            f.write( u'\t%s'%pathway_sel[(name_pathway,sga)])
        f.write(u'\n')
    f.close()

    # f.write( u'%s(%s,Y)\n'%(relname,src) )

    # f = open(path_outputData+'/node.txt', 'w')
    # print >> f, 'Id\tTP53pathway'
    # for line in tp53_sel.keys():
    #     print >> f, line + '\t'+ str(tp53_sel[line])
    # f.close()


    f = open(path_outputData+'/linked.txt', 'w')
    print >> f, 'Source\tTarget\tType\tWeight'
    for line in sgaWsga.keys():
        sga1, sga2 = line[0], line[1]
        print >> f, sga1+'\t'+sga2+'\tUndirected\t'+str(sgaWsga[(sga1,sga2)])
    f.close()


#EOF.
