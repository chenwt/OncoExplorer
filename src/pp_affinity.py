# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd
import numpy as np
import scipy.io
import os

# def getAffinity(cancer):
#     set_graph = set()
#     tmp_graph = dd(float)
#     #print 'reading from... {}'.format('/usr1/public/yifeng/Github/inputData/ensemble.txt')
#     k = 0
#     flag = True
#     sga2pat = dd(set) # Tumors that the sga exists.
#     #sga2patnum = dd()
#     with open('/usr1/public/yifeng/Github/inputData/ensemble.txt') as f:
#         next(f)
#         for line in f:
#             line = line.strip().split('\t')
#             patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
#             if cancer == 'pancan': # pancan is alway true
#                 flag = True
#             elif can == cancer: # cancer name match specific can
#                 flag = True
#             else: # cancer does not match the specific can
#                 flag = False
#             if not flag: continue
#             k = k+1
#             #print patid, can, sga, deg, prob
#             #if k%100000 == 0: print k 
#             if (sga,deg) not in set_graph:
#                 tmp_graph[(sga,deg)] = 1.0
#                 set_graph.add((sga,deg))
#             else:
#                 tmp_graph[(sga,deg)] += 1.0
#             sga2pat[sga].add(patid)
#     print '#readin_edge of {} = {}'.format(cancer, k)

#     for sga in sga2pat.keys():
#         print sga+'\t'+str(len(sga2pat[sga]))

#     f = open('/usr1/public/yifeng/Github/OncoExplorer/src/bond_'+cancer+'.txt', 'w')
#     for line in tmp_graph.keys():
#         sga,deg = line[0], line[1]
#         print >> f, str(tmp_graph[line])
#     f.close()

#     # Calculate the overlap of DEGs between SGAs.
#     #sga2deg = dd(list)
#     # There should be a threshold here to filter the set_graph.
def getPathway(pathway):
    set_pathway = set()
    path = '/usr1/public/yifeng/Github/inputData/kegg_html/'+pathway+'_html.txt'
    f = open(path, 'r')
    for line in f:
        line = line.strip().split('\t')
        genes = line[3]
        genes = genes.split('\"')[1].split(', ')
        for gene in genes:
            gene = gene.split(' (')[1].split(')')[0].lower()
            set_pathway.add(gene)
    f.close()
    return set_pathway

def getPathwayGo(pathway):
    set_pathway = set()
    path = '/usr1/public/yifeng/Github/inputData/kegg_html/'+pathway+'_go.txt'

    if not os.path.exists(path): return set_pathway

    f = open(path, 'r')
    for line in f:
        gene = line.strip().split('\t')[1].lower()
        gene = gene.split('(')[0]
        #print gene
        set_pathway.add(gene)
    f.close()
    return set_pathway

def getCoordinates(cancer, set_notch1, set_p53, set_pi3k):
    # sga -> set of degs
    sga2deg = dd(set)
    tmp_graph = dd(float)
    # set of DEGs that appear in the (SGA, DEG) pairs
    set_deg = set()
    #print 'reading from... {}'.format('/usr1/public/yifeng/Github/inputData/ensemble.txt')
    k = 0
    flag = True
    sga2pat = dd(set) # Tumors that the sga exists.
    #sga2patnum = dd()
    with open('/usr1/public/yifeng/Github/inputData/ensemble.txt') as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            patid, can, sga, deg, prob = int(line[1]),line[2],line[4],line[5],line[6]
            if cancer == 'pancan': # pancan is alway true
                flag = True
            elif can == cancer: # cancer name match specific can
                flag = True
            else: # cancer does not match the specific can
                flag = False
            if not flag: continue
            k = k+1
            #print patid, can, sga, deg, prob
            set_deg.add(deg)
            if deg not in sga2deg[sga]:
                tmp_graph[(sga,deg)] = 1.0
                sga2deg[sga].add(deg)
            else:
                tmp_graph[(sga,deg)] += 1.0
            sga2pat[sga].add(patid)
    print '#readin_edge of {} = {}'.format(cancer, k)

    fi = open('/usr1/public/yifeng/Github/OncoExplorer/src/index2name_'+cancer+'.txt', 'w')
    print >> fi, 'index\tgene\tnotch1\tp53\tpi3k'

    set_deg = list(set_deg)
    set_sga = list(sga2deg.keys())
    corr = np.zeros([len(set_sga), len(set_deg)], float)

    iter_sga = 0

    for sga in set_sga:
        iter_sga += 1
        innotch1 = 0
        inpi3k = 0
        inp53 = 0
        if sga in set_notch1: innotch1 = 1
        if sga in set_p53: inp53 = 1
        if sga in set_pi3k: inpi3k = 1
        print >> fi, str(iter_sga)+'\t'+sga+'\t'+str(innotch1)+'\t'+str(inp53)+'\t'+str(inpi3k)
        numOfPat = len(sga2pat[sga])
        iter_deg = 0
        for deg in set_deg:
            iter_deg += 1
            #if (sga, deg) in tmp_graph.keys():
            corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]/numOfPat

    fi.close()
    scipy.io.savemat('/usr1/public/yifeng/Github/OncoExplorer/src/bond_'+cancer+'.mat', mdict = {'corr': corr})


if __name__ == '__main__':


            #print gene
        #print genes
    set_notch1 = getPathwayGo('notch1')
    set_p53 = getPathwayGo('p53')
    set_pi3k = getPathwayGo('pi3k')

    print len(set_notch1), len(set_p53), len(set_pi3k)

    for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        getCoordinates(cancer, set_notch1, set_p53, set_pi3k)

#EOF.
