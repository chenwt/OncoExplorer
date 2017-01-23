# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd
import numpy as np
import scipy.io

def getAffinity(cancer):
    set_graph = set()
    tmp_graph = dd(float)
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
            #if k%100000 == 0: print k 
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            else:
                tmp_graph[(sga,deg)] += 1.0
            sga2pat[sga].add(patid)
    print '#readin_edge of {} = {}'.format(cancer, k)

    for sga in sga2pat.keys():
        print sga+'\t'+str(len(sga2pat[sga]))

    f = open('/usr1/public/yifeng/Github/OncoExplorer/src/bond_'+cancer+'.txt', 'w')
    for line in tmp_graph.keys():
        sga,deg = line[0], line[1]
        print >> f, str(tmp_graph[line])
    f.close()

    # Calculate the overlap of DEGs between SGAs.
    #sga2deg = dd(list)
    # There should be a threshold here to filter the set_graph.

def getCoordinates(cancer):
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
            #if k%100000 == 0: print k 
            set_deg.add(deg)
            if deg not in sga2deg[sga]:
                tmp_graph[(sga,deg)] = 1.0
                sga2deg[sga].add(deg)
            else:
                tmp_graph[(sga,deg)] += 1.0
            sga2pat[sga].add(patid)
    print '#readin_edge of {} = {}'.format(cancer, k)

    # for sga in sga2pat.keys():
    #     print sga+'\t'+str(len(sga2pat[sga]))

    #f = open('/usr1/public/yifeng/Github/OncoExplorer/src/bond_'+cancer+'.txt', 'w')
    fi = open('/usr1/public/yifeng/Github/OncoExplorer/src/index2name_'+cancer+'.txt', 'w')
    print >> fi, 'index\tgene'

    set_deg = list(set_deg)
    set_sga = list(sga2deg.keys())
    #corr = [[0 for y in range(len(set_deg))] for x in range(len(sga2deg.keys()))] 
    corr = np.zeros([len(set_sga), len(set_deg)], float)

    iter_sga = 0

    for sga in set_sga:
        iter_sga += 1
        print >> fi, str(iter_sga)+'\t'+sga
        numOfPat = len(sga2pat[sga])
        iter_deg = 0
        #print iter_sga
        for deg in set_deg:
            iter_deg += 1
            #if (sga, deg) in tmp_graph.keys():
            corr[iter_sga-1][iter_deg-1] = 1.0*tmp_graph[(sga,deg)]/numOfPat
            #else:
            #    corr[iter_sga-1][iter_deg-1] = 0.0

    #for x in range(len(sga2deg.keys())):
        #for y in range(len(set_deg)):
    #        print >> f, corr[x]

    #set_sga = 
    # for line in tmp_graph.keys():
    #     sga,deg = line[0], line[1]
    #     print >> f, str(tmp_graph[line])
    #f.close()
    fi.close()
    scipy.io.savemat('/usr1/public/yifeng/Github/OncoExplorer/src/bond_'+cancer+'.mat', mdict = {'corr': corr})

    #return corr

if __name__ == '__main__':

    # Getting the isA.cfacts file.
    # path_inputData = '/usr1/public/yifeng/Github/inputData'

    # path_outputData = '/usr1/public/yifeng/Github/outputData'

    for cancer in ['ov']:
    #for cancer in ['pancan', 'brca', 'gbm', 'ov']:
        getCoordinates(cancer)

    # path_outputData = '/usr1/public/yifeng/Github/PropprData'
    # f = open(path_outputData+'/Drives_'+cancer+'.graph', 'w')
    # print >> f, 'SGA\tDEG\tWeight'
    # for line in tmp_graph.keys():
    #     sga,deg,weight = line[0],line[1],tmp_graph[(line[0],line[1])]
    #     print >> f, sga+'\t'+deg+'\t'+str(weight)
    #print >> f, 'SGA\tDEG'



    # for line in tmp_graph.keys():
    #     sga,deg,weight = line[0],line[1],tmp_graph[(line[0],line[1])]
    #     print sga+'\t'+deg+'\t'+str(weight)



    #f.close()
#EOF.
