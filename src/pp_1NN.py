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
        sga2sgaAaff = getAffinityBaseline1(cancer)

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
