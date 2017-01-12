from collections import defaultdict as dd
import io
import os
import random
import argparse


if __name__ == '__main__':
    ''' 
    Group the SGAs and generate files for matlab.
    '''
    cancer = 'ov'
    path_sol = '/usr1/public/yifeng/Github/PropprData/test_'+cancer+'.solutions.txt'

    gene_set = set()
    f = open(path_sol, 'r')
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            line = line.split('\t')[1].split('(')[1].split(')')[0].split(',')[0]
            gene_set.add(line)

    f.close()



    gene2base = dd()
    for line in gene_set:
        gene2base[line] = 'undifined'

    f = open(path_sol, 'r')
    for line in f:
        line = line.strip()
        if line.startswith('1\t'):
            gene = line.split('\t')[2].split(',')[0].split('(')[1]
            base = line.split('\t')[2].split(',')[1].split(')')[0][4:]
            #line = line.split('\t')[1].split('(')[1].split(')')[0].split(',')[0]
            gene2base[gene] = base
            #print gene, base
            #gene_set.add(line)

    f.close()


    # path_sol = '/usr1/public/yifeng/Github/outputData/test_'+cancer+'.txt'
    # f = open(path_sol, 'w')
    # print >> f, 'gene\tgoterm'
    # for line in gene_set:
    #     print >> f, line+'\t'+gene2base[line]

    # f.close()

    base_set = set()
    base_set.add('undifined')
    path_base = '/usr1/public/yifeng/Github/PropprData/isBase.cfacts'
    f = open(path_base, 'r')
    for line in f:
        line = line.strip().split('\t')[1]
        base_set.add(line)

    f.close()

    goterm2group = dd()

    #path_lut = '/usr1/public/yifeng/Github/outputData/lut.txt'

    #f = open(path_lut, 'w')
    #print >> f, 'Id\tGroup'
    k = 0
    for line in base_set:
        k += 1
        goterm2group[line] = k

    f.close()


    path_sol = '/usr1/public/yifeng/Github/outputData/test_'+cancer+'.txt'
    f = open(path_sol, 'w')
    print >> f, 'gene\tgoterm\tgroup'
    for line in gene_set:
        print >> f, line+'\t'+gene2base[line]+'\t'+str(goterm2group[gene2base[line]])

    f.close()

    print 'Done!'

#EOF.
