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


def test(gene2goId, set_pathway, choice, metric):

    if choice == 'Proppr':

        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            sga2sgaAaff = getAffinityProppr(cancer)

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
                    # if ovlp > 0:
                    #     count_true += 1
                        #print sga, NNgene, NNdist
                acc = 1.0*count_true/count_total
                avgacc += acc
                # Omit if necessary
                #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
            print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)
    #getAffinityProppr


    if choice == 'Lu affinity':
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            sga2sgaAaff = getLuAffinity(cancer)

            #for pathway in ['notch1','p53','pi3k']:
            #set_pathway = getPathway(pathway)

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


    if choice == 'Lu dist':
        print metric, choice
        for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        #for cancer in ['brca']:
            sga2sgaDist = getCoordinates(cancer)

            #for pathway in ['notch1','p53','pi3k']:
            #    set_pathway = getPathway(pathway)

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
        #for cancer in ['brca']:
            sga2sgaAaff = getAffinityBaseline1(cancer)

            #for pathway in ['notch1','p53','pi3k']:
                #set_pathway = getPathway(pathway)

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


    # choice = 'Proppr'
    #choice = 'hello world'

    gene2goId, set_pathway = getGene2GoId()

    #print len(gene2goId), len(set_pathway)

    path = '../lab/gene2goId.txt'
    f = open(path, 'w')
    for g in gene2goId.keys():
        for i in gene2goId[g]:
            print >> f, g+'\t'+i
    f.close()


    #print len(gene2goId), len(set_pathway)

    # metric = 'one'
    # metric = 'total'
    # metric = 'jaccard'

#    for metric in ['one', 'total', 'jaccard']:
#        for choice in ['Random', 'Baseline1', 'Proppr', 'Lu affinity', 'Jaccard1', 'Lu dist']:
#            test(gene2goId, set_pathway, choice, metric)
    metric = 'one'
    choice = 'Baseline1'

    print metric, choice
    for cancer in ['pancan', 'brca', 'gbm', 'ov']:
    #for cancer in ['brca']:
        sga2sgaAaff = getAffinityBaseline1(cancer)

        path = '../lab/sga2sgaAff_baseline_'+cancer+'.txt'
        f = open(path, 'w')
        print >> f, 'Source\tTarget\tWeight'
        for sga1 in sga2sgaAaff.keys():
            for line in sga2sgaAaff[sga1]:
                sga2, overlap = line[0], line[1]
                if overlap == 0: continue
                if sga1 == sga2: continue
                print >> f, sga1 + '\t' + sga2 + '\t' + str(overlap)
                #print >> f, sga1 + '\t' + sga2 + '\t' + overlap
        f.close()

#sga2sgaAaff[sga1].add((sga2, overlap))


        #for pathway in ['notch1','p53','pi3k']:
        #    set_pathway = getPathway(pathway)

        # avgacc = 0
        # for trial in range(10):
        #     count_total = 0
        #     count_true = 0
        #     for sga in set_pathway:
        #         if sga not in sga2sgaAaff.keys(): continue
        #         count_total += 1
        #         NNgene = 'helloworld'
        #         NNdist = 0.0
        #         Context = list(sga2sgaAaff[sga])
        #         shuffle(Context)
        #         for line in Context:
        #             sgaContext, aff = line[0], line[1]
        #             if sgaContext == sga: continue
        #             if aff > NNdist:
        #                 NNdist = aff
        #                 NNgene = sgaContext
        #         #print sga, NNgene, NNdist

        #         ovlp = gene2goId[sga].intersection(gene2goId[NNgene])
        #         ovlp = len(ovlp)

        #         if metric == 'one':
        #             if ovlp > 0:
        #                 count_true += 1
        #         elif metric == 'total':
        #             count_true += ovlp
        #         elif metric == 'jaccard':
        #             count_true += 1.0*ovlp/len(gene2goId[sga].union(gene2goId[NNgene]))
        #     acc = 1.0*count_true/count_total
        #     avgacc += acc
        #     # Omit if necessary
        #     #print '\tcancer:'+cancer+'\ttrial:'+str(trial)+'\taccuracy:'+str(acc)
        # print 'cancer:'+cancer+'\tavgacc:'+str(1.0*avgacc/10)+'\t\t#checked_genes:'+str(count_total)

#EOF.
