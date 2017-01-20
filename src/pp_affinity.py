# Get the affinity between SGAs of specific tumor.
from collections import defaultdict as dd

def getAffinity(cancer):
    set_graph = set()
    tmp_graph = dd()
    #print 'reading from... {}'.format('/usr1/public/yifeng/Github/inputData/ensemble.txt')
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
            #print patid, can, sga, deg, prob
            #if k%100000 == 0: print k 
            if (sga,deg) not in set_graph:
                tmp_graph[(sga,deg)] = 1.0
                set_graph.add((sga,deg))
            else:
                tmp_graph[(sga,deg)] += 1.0
    print '#readin_edge of {} = {}'.format(cancer, k)

if __name__ == '__main__':

    # Getting the isA.cfacts file.
    # path_inputData = '/usr1/public/yifeng/Github/inputData'

    # path_outputData = '/usr1/public/yifeng/Github/outputData'

    cancer = 'brca'
    getAffinity(cancer)
    cancer = 'gbm'
    getAffinity(cancer)
    cancer = 'ov'
    getAffinity(cancer)
    cancer = 'pancan'
    getAffinity(cancer)
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
