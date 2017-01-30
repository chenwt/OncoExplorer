# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd
import numpy as np
import scipy.io
import os
from random import shuffle


def getAffinityBaseline1(cancer):
    sga2sgaAaff = dd(set)
    sga2deg = dd(set)
    k = 0
    flag = True
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

            sga2deg[sga].add(deg)

    #print '#readin_edge of {} = {}'.format(cancer, k)

    f = open('/usr1/public/yifeng/Github/anal/corr_'+cancer+'_baseline1.txt', 'w')
    #print >> f, 'Source\tTarget\tWeight'
    #traversedset = set()
    k = 0
    for sga1 in sga2deg.keys():
        #traversedset.add(sga1)
        for sga2 in sga2deg.keys():
            #if sga2 != sga1:
            #    if sga2 in traversedset: continue 
            
            overlap = len(sga2deg[sga1].intersection(sga2deg[sga2]))
            #if overlap == 0: continue
            k += 1
            print >> f, sga1+'\t'+sga2+'\t'+str(overlap)
            sga2sgaAaff[sga1].add((sga2, overlap))
        #degset = sga2deg[sga]
    numsga = len(sga2deg.keys())

    #print numsga, numsga*numsga, k
    #print len(sga2sgaAaff), len(sga2sgaAaff.keys())

    f.close()
    return sga2sgaAaff


def getAffinityJaccard1(cancer):
    sga2sgaAaff = dd(set)
    sga2deg = dd(set)
    k = 0
    flag = True
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
            sga2deg[sga].add(deg)
    #print '#readin_edge of {} = {}'.format(cancer, k)


    f = open('/usr1/public/yifeng/Github/anal/corr_'+cancer+'_jaccard1.txt', 'w')
    print >> f, 'Source\tTarget\tWeight'
    #traversedset = set()
    k = 0
    for sga1 in sga2deg.keys():
        #traversedset.add(sga1)
        for sga2 in sga2deg.keys():
            #if sga2 != sga1:
            #    if sga2 in traversedset: continue 
            
            overlap = 1.0*len(sga2deg[sga1].intersection(sga2deg[sga2]))/len(sga2deg[sga1].union(sga2deg[sga2]))
            #if overlap < 1.0*threshold/100: continue
            k += 1
            print >> f, sga1+'\t'+sga2+'\t'+str(overlap)
            sga2sgaAaff[sga1].add((sga2, overlap))
        #degset = sga2deg[sga]
    #numsga = len(sga2deg.keys())

    #print numsga, numsga*numsga, k
    # for line in tmp_graph.keys():
    #     sga,deg = line[0], line[1]
    #     print >> f, str(tmp_graph[line])
    f.close()

    return sga2sgaAaff


def getCoordinatesCount(cancer, set_notch1, set_p53, set_pi3k):
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
            corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]

    fi.close()
    scipy.io.savemat('/usr1/public/yifeng/Github/OncoExplorer/src/bond_'+cancer+'.mat', mdict = {'corr': corr})

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


if __name__ == '__main__':


    #print len(set_notch1), len(set_p53), len(set_pi3k)


    for cancer in ['pancan', 'brca', 'gbm', 'ov']:
    #for cancer in ['brca']:
        sga2sgaAaff = getAffinityJaccard1(cancer)

        for pathway in ['notch1','p53','pi3k']:
            set_pathway = getPathway(pathway)

            avgacc = 0
            for trial in range(10):
                count_total = 0
                count_true = 0
                for sga in set_pathway:
                    if sga not in sga2sgaAaff.keys(): continue
                    count_total += 1
                    NNgene = 'helloworld'
                    NNdist = 0.0
                    Context = list(sga2sgaAaff[sga])
                    shuffle(Context)
                    for line in Context:
                        sgaContext, aff = line[0], line[1]
                        if sgaContext == sga: continue
                        if aff > NNdist:
                            NNdist = aff
                            NNgene = sgaContext
                    #print sga, NNgene, NNdist
                    if NNgene in set_pathway:
                        count_true += 1
                        #print sga, NNgene, NNdist
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print pathway+'\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print pathway+'\tcancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)
                #print sga, sgaContext, aff




        #for threshold in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        #    getAffinityJaccard1(cancer, threshold)



    #for cancer in ['pancan', 'brca', 'gbm', 'ov']:
    #    getCoordinatesCount(cancer, set_notch1, set_p53, set_pi3k)

#EOF.
