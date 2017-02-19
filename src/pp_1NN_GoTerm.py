# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd
import numpy as np
import scipy.io
import os
from random import shuffle
import math

# Get pathway information
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


# Return: Gene->GoTerm, genes that have BP GoTerm
def getGene2GoId():
    goAnn = set()
    path_goAnn = '/usr1/public/yifeng/Github/inputData/go/goa_human.gaf'
    f = open(path_goAnn, 'r')
    for line in f:
        if line[0] == '!': continue
        line = line.split('\t')
        gene, goId = line[2], line[4]
        gene = gene.lower()
        goId = 'go'+goId.split(':')[1].lower()
        goAnn.add((gene, goId))
    f.close()

    f = open('/usr1/public/yifeng/Github/inputData/go/go-basic.obo', 'r')
    go_bp_set = set() # The set that belongs to biological process.

    flag_term = False
    src_id = ''
    name_space = ''
    for line in f:
        line = line.strip()
        if line == '[Term]': flag_term = True

        if flag_term == False: continue

        if line.startswith('id: GO:'):
            src_id = line
            continue
        if line.startswith('namespace:'):
            name_space = line
            continue
        if line == '':
            if name_space == 'namespace: biological_process':
                src = 'go'+src_id.split(':')[2]
                go_bp_set.add(src)

            flag_term = False
            continue
    f.close()

    gene2goId = dd(set)
    set_gene = set()
    for line in goAnn:
        gene, goId = line[0], line[1]
        if goId in go_bp_set:
            gene2goId[gene].add(goId)
            set_gene.add(gene)

    return gene2goId, set_gene



# Get information of affinity
# return: sga2sgaAaff sga->(sga, aff)

# version of max
def getAffinityProppr(cancer):
    sga2sgaAaff = dd(set)
    sgaAsga2prob = dd(float)

    set_sga = set()
    f = open('/usr1/public/yifeng/Github/PropprData/'+cancer+'.solutions.txt', 'r')
    for line in f:
        if line.startswith('#'): continue
        line1 = line.strip().split('\t')[-1].split('(')[1].split(')')[0].split(',')
        prob = line.strip().split('\t')[1]
        sga1, sga2, prob = line1[0], line1[1], float(prob)
        set_sga.add(sga1)
        set_sga.add(sga2)
        sgaAsga2prob[(sga1,sga2)] = prob
    f.close()

    for sga1 in set_sga:
        for sga2 in set_sga:
            #prob = sgaAsga2prob[(sga1,sga2)]*sgaAsga2prob[(sga2,sga1)]
            #prob = sgaAsga2prob[(sga1,sga2)]+sgaAsga2prob[(sga2,sga1)]
            prob = max(sgaAsga2prob[(sga1,sga2)], sgaAsga2prob[(sga2,sga1)])
            
            sga2sgaAaff[sga1].add((sga2,prob))
    return sga2sgaAaff


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

def getAffinityDW(wl,nw):
    #sga2sgaAaff = dd(set)
    #sga2deg = dd(set)

    # gid -> gene
    gid2gene = dd()
    path = '/usr1/public/yifeng/Github/deepwalk-master/data/LUT_DW_brca.txt'
    f = open(path, 'r')
    rg_gid = 0
    for line in f:
        
        l = line.strip().split('\t')
        gene, gid = l[0], int(l[1])
        
        if gene.endswith('_deg'): break
        rg_gid += 1
        gid2gene[rg_gid-1] = gene
        #print gene, gid, rg_gid
    f.close()

    W = np.zeros([rg_gid, 40], float)
    path = '/usr1/public/yifeng/Github/deepwalk-master/result/brca_km40_wl'+str(wl)+'_nw'+str(nw)+'.embeddings'
    f = open(path, 'r')
    next(f)
    k = 0
    for line in f:
        l = line.strip().split(' ', 1)
        gid = int(l[0])
        if gid <= rg_gid:
            k += 1
            x = [float(i) for i in l[1].split(' ')]
            W[gid-1,:] = x
            #print k, gid, x
    f.close()
    
    sga2sgaAaff = dd(set)
    
    M = 0
    for id1 in range(rg_gid):
        g1 = gid2gene[id1]
        for id2 in range(rg_gid):
            g2 = gid2gene[id2]
            d = W[id1,:] - W[id2,:]
            d =  np.linalg.norm(d)
            M = max(d, M)
            
    for id1 in range(rg_gid):
        g1 = gid2gene[id1]
        for id2 in range(rg_gid):
            g2 = gid2gene[id2]
            d = W[id1,:] - W[id2,:]
            d =  M - np.linalg.norm(d)
            
            sga2sgaAaff[g1].add((g2, d))
    
    return sga2sgaAaff

def getAffinityWeightedOverlap(cancer, Mmpa):
    sga2sgaAaff = dd(set)
    sgaAdeg2w = dd(float)
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
            sgaAdeg2w[(sga,deg)] += 1.0
            sga2deg[sga].add(deg)

    #print '#readin_edge of {} = {}'.format(cancer, k)

    k = 0
    for sga1 in sga2deg.keys():
        for sga2 in sga2deg.keys():
            set_CommonDeg = sga2deg[sga1].intersection(sga2deg[sga2])
            overlap = 0.0
            for d in set_CommonDeg:
                if Mmpa == 'prod': overlap += sgaAdeg2w[(sga1,d)]*sgaAdeg2w[(sga2,d)]
                if Mmpa == 'min': overlap += min(sgaAdeg2w[(sga1,d)],sgaAdeg2w[(sga2,d)])
                if Mmpa == 'max': overlap += max(sgaAdeg2w[(sga1,d)],sgaAdeg2w[(sga2,d)])
                if Mmpa == 'avg': overlap += 0.5*(sgaAdeg2w[(sga1,d)]+sgaAdeg2w[(sga2,d)])
            #overlap = len()
            #if overlap == 0: continue
            k += 1
            #print sga1+'\t'+sga2+'\t'+str(overlap)
            sga2sgaAaff[sga1].add((sga2, overlap))
        #degset = sga2deg[sga]
    #numsga = len(sga2deg.keys())

    #print numsga, numsga*numsga, k
    #print len(sga2sgaAaff), len(sga2sgaAaff.keys())


    return sga2sgaAaff

#TODO:
def getAffinityTFIDF(cancer, Mmpa):
    sga2sgaAaff = dd(set)
    deg2sga = dd(set)
    degSpan = dd(float)
    
    sga2deglist = dd(list)
    #sgaAdeg2w = dd(float)
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
            deg2sga[deg].add(sga)
            sga2deglist[sga].append(deg)
            #sgaAdeg2w[(sga,deg)] += 1.0
            sga2deg[sga].add(deg)


    Nsga = len(sga2deg.keys())
    #print Nsga
    for deg in deg2sga.keys():
        degSpan[deg] = len(deg2sga[deg])
        #print 'degSpan',deg, degSpan[deg]
    #print '#readin_edge of {} = {}'.format(cancer, k)

    k = 0
    for sga1 in sga2deg.keys():
        for sga2 in sga2deg.keys():
            set_CommonDeg = sga2deg[sga1].intersection(sga2deg[sga2])
            overlap = 0.0
            #print sga1, sga2, overlap
            for d in set_CommonDeg:
                
                sga1d_tf = 1.0*sga2deglist[sga1].count(d)/len(sga2deglist[sga1])
                sga1d_idf = 1.0*math.log(1.0*Nsga/degSpan[d])
                
                sga1d_tfidf = 1.0*sga1d_tf*sga1d_idf
                sga2d_tf = 1.0*sga2deglist[sga2].count(d)/len(sga2deglist[sga2])
                sga2d_idf = 1.0*math.log(1.0*Nsga/degSpan[d])
                sga2d_tfidf = 1.0*sga2d_tf*sga2d_idf
                #print sga2deglist[sga1]
                #print d
                #print sga1d_tf, sga2deglist[sga1].count(d), sga1d_idf, sga2d_tf, sga2deglist[sga2].count(d), sga2d_idf
                if Mmpa == 'prod': overlap += sga1d_tfidf*sga2d_tfidf
                if Mmpa == 'min': overlap += min(sga1d_tfidf, sga2d_tfidf)
                if Mmpa == 'max': overlap += max(sga1d_tfidf, sga2d_tfidf)
                if Mmpa == 'avg': overlap += 0.5*(sga1d_tfidf+sga2d_tfidf)
                #print sga1, sga2, overlap
            #overlap = len()
            #if overlap == 0: continue
            k += 1
            #print sga1+'\t'+sga2+'\t'+str(overlap)
            
            sga2sgaAaff[sga1].add((sga2, overlap))
        #degset = sga2deg[sga]
    #numsga = len(sga2deg.keys())

    #print numsga, numsga*numsga, k
    #print len(sga2sgaAaff), len(sga2sgaAaff.keys())


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
            #print >> f, sga1+'\t'+sga2+'\t'+str(overlap)
            sga2sgaAaff[sga1].add((sga2, overlap))
        #degset = sga2deg[sga]
    #numsga = len(sga2deg.keys())

    #print numsga, numsga*numsga, k
    # for line in tmp_graph.keys():
    #     sga,deg = line[0], line[1]
    #     print >> f, str(tmp_graph[line])
    f.close()

    return sga2sgaAaff

def getCoordinates(cancer):

    sga2sgaAaff = dd(set)
    sga2deg = dd(set)

    tmp_graph = dd(float)
    # set of DEGs that appear in the (SGA, DEG) pairs
    set_deg = set()
    k = 0
    flag = True
    sga2pat = dd(set) # Tumors that the sga exists.
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
    #print '#readin_edge of {} = {}'.format(cancer, k)


    set_deg = list(set_deg)
    set_sga = list(sga2deg.keys())
    corr = np.zeros([len(set_sga), len(set_deg)], float)

    #corr = dd(list)
    iter_sga = 0

    for sga in set_sga:
        iter_sga += 1
        #if iter_sga % 10 == 0: print iter_sga
        numOfPat = len(sga2pat[sga])
        iter_deg = 0
        for deg in set_deg:
            iter_deg += 1
            #if (sga, deg) in tmp_graph.keys():
            corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]/numOfPat
                #corr[sga].append(1.0*tmp_graph[(sga,deg)]/numOfPat)

    iter_sga1 = 0
    for sga1 in set_sga:
        iter_sga1 += 1
        #if iter_sga1 % 10 == 0: print iter_sga1
        #traversedset.add(sga1)
        v1 = corr[iter_sga1-1]
        iter_sga2 = 0
        for sga2 in set_sga:
            iter_sga2 += 1
            v2 = corr[iter_sga2-1]
            #if sga2 != sga1:
            #    if sga2 in traversedset: continue 
            
            #aff = - distSq(corr[sga1], corr[sga2])
            aff1 = v1 - v2
            aff = np.inner(aff1, aff1)
            #overlap = 1.0*len(sga2deg[sga1].intersection(sga2deg[sga2]))/len(sga2deg[sga1].union(sga2deg[sga2]))
            #if overlap < 1.0*threshold/100: continue
            #k += 1
            #print >> f, sga1+'\t'+sga2+'\t'+str(overlap)
            sga2sgaAaff[sga1].add((sga2, aff))
        #print sga, len(corr[sga])
            #corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]/numOfPat
    return sga2sgaAaff


def getLuAffinity(cancer):

    sga2sgaAaff = dd(set)
    sga2deg = dd(set)

    tmp_graph = dd(float)
    # set of DEGs that appear in the (SGA, DEG) pairs
    set_deg = set()
    k = 0
    flag = True
    sga2pat = dd(set) # Tumors that the sga exists.
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
    #print '#readin_edge of {} = {}'.format(cancer, k)


    set_deg = list(set_deg)
    set_sga = list(sga2deg.keys())
    corr = np.zeros([len(set_sga), len(set_deg)], float)

    #corr = dd(list)
    iter_sga = 0

    for sga in set_sga:
        iter_sga += 1
        #if iter_sga % 10 == 0: print iter_sga
        numOfPat = len(sga2pat[sga])
        iter_deg = 0
        for deg in set_deg:
            iter_deg += 1
            #if (sga, deg) in tmp_graph.keys():
            corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]/numOfPat
                #corr[sga].append(1.0*tmp_graph[(sga,deg)]/numOfPat)

    iter_sga1 = 0
    for sga1 in set_sga:
        iter_sga1 += 1
        #if iter_sga1 % 10 == 0: print iter_sga1
        #traversedset.add(sga1)
        v1 = corr[iter_sga1-1]
        iter_sga2 = 0
        for sga2 in set_sga:
            iter_sga2 += 1
            v2 = corr[iter_sga2-1]
            #if sga2 != sga1:
            #    if sga2 in traversedset: continue 
            
            #aff = - distSq(corr[sga1], corr[sga2])
            #aff1 = v1 - v2
            aff = np.inner(v1, v2)
            #overlap = 1.0*len(sga2deg[sga1].intersection(sga2deg[sga2]))/len(sga2deg[sga1].union(sga2deg[sga2]))
            #if overlap < 1.0*threshold/100: continue
            #k += 1
            #print >> f, sga1+'\t'+sga2+'\t'+str(overlap)
            sga2sgaAaff[sga1].add((sga2, aff))
        #print sga, len(corr[sga])
            #corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]/numOfPat
    return sga2sgaAaff










def test(gene2goId, set_pathway, choice, metric):

    if choice == 'Proppr':

        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
            sga2sgaAaff = getAffinityProppr(cancer)
            
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
                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)
                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)

    if choice == 'Lu affinity':
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
            sga2sgaAaff = getLuAffinity(cancer)
            avgacc = 0
            for trial in range(10):
                count_total = 0
                count_true = 0
                for sga in set_pathway:
                    if sga not in sga2sgaAaff.keys(): continue
                    count_total += 1
                    NNgene = 'helloworld'
                    NNdist = -10000.0
                    Context = list(sga2sgaAaff[sga])
                    shuffle(Context)
                    for line in Context:
                        sgaContext, aff = line[0], line[1]
                        if sgaContext == sga: continue
                        if aff > NNdist:
                            NNdist = aff
                            NNgene = sgaContext

                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                        #print sga, NNgene, NNdist
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)


    if choice == 'Lu dist':
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
            sga2sgaDist = getCoordinates(cancer)

            avgacc = 0
            for trial in range(10):
                count_total = 0
                count_true = 0
                for sga in set_pathway:
                    if sga not in sga2sgaDist.keys(): continue
                    count_total += 1
                    NNgene = 'helloworld'
                    NNdist = 1000000000000.0
                    Context = list(sga2sgaDist[sga])
                    shuffle(Context)
                    for line in Context:
                        sgaContext, dist = line[0], line[1]
                        if sgaContext == sga: continue
                        if dist < NNdist:
                            NNdist = dist
                            NNgene = sgaContext
                    #print sga, NNgene, NNdist


                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                        #print sga, NNgene, NNdist
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print 'cancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)
                #print sga, sgaContext, aff


    if choice == 'Random':
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
            sga2sgaAaff = getAffinityBaseline1(cancer)

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
                        #if aff > NNdist:
                        NNdist = aff
                        NNgene = sgaContext
                    #print sga, NNgene, NNdist
                    
                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                        #print sga, NNgene, NNdist
                acc = 1.0*count_true/count_total
                avgacc += acc
                    # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)
    if choice == 'Baseline1':
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            sga2sgaAaff = getAffinityBaseline1(cancer)

            #for pathway in ['notch1','p53','pi3k']:
            #    set_pathway = getPathway(pathway)

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

                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)


#TODO:DW
    if choice == 'DW':
        print metric, choice
        for cancer in ['brca']:
        #for cancer in ['brca']:
            for wl in [1,2,3,5,10,20,40]:
                for nw in [10,20,50,100,200]:

                    sga2sgaAaff = getAffinityDW(wl,nw)
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
        
                            ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                            ovlp = len(ovlp)
        
                            if metric == 'one':
                                if ovlp > 0:
                                    count_true += 1
                            elif metric == 'total':
                                count_true += ovlp
                            elif metric == 'jaccard':
                                count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                        acc = 1.0*count_true/count_total
                        avgacc += acc
                        # Omit if necessary
                        #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
                    print 'cancer:'+cancer+'\twl:'+str(wl)+'\tnw:'+str(nw)+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)

    if choice in ['WeightedOverlap_min', 'WeightedOverlap_max', 'WeightedOverlap_prod','WeightedOverlap_avg']:
        if choice == 'WeightedOverlap_min': Mmpa = 'min'
        if choice == 'WeightedOverlap_max': Mmpa = 'max'
        if choice == 'WeightedOverlap_prod': Mmpa = 'prod'
        if choice == 'WeightedOverlap_avg': Mmpa = 'avg'
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            
            #Mmpa = 'min'
            sga2sgaAaff = getAffinityWeightedOverlap(cancer, Mmpa)

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

                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)

    #TODO:
    if choice in ['TFIDF_min','TFIDF_max','TFIDF_avg','TFIDF_prod']:
        if choice == 'TFIDF_min': Mmpa = 'min'
        if choice == 'TFIDF_max': Mmpa = 'max'
        if choice == 'TFIDF_prod': Mmpa = 'prod'
        if choice == 'TFIDF_avg': Mmpa = 'avg'
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            
            #Mmpa = 'min'
            sga2sgaAaff = getAffinityTFIDF(cancer, Mmpa)

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

                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)


    if choice == 'Jaccard1':

        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            sga2sgaAaff = getAffinityJaccard1(cancer)

            #for pathway in ['notch1','p53','pi3k']:
            #    set_pathway = getPathway(pathway)

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
                    ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
                    ovlp = len(ovlp)

                    if metric == 'one':
                        if ovlp > 0:
                            count_true += 1
                    elif metric == 'total':
                        count_true += ovlp
                    elif metric == 'jaccard':
                        count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))

                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)

if __name__ == '__main__':

    gene2goId, set_pathway = getGene2GoId()

    for metric in ['one', 'total', 'jaccard']:
        for choice in ['DW']:
#    for metric in ['one']:
#        for choice in ['TFIDF_min']:            
            test(gene2goId, set_pathway, choice, metric)



#EOF.
