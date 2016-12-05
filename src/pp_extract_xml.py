from collections import defaultdict as dd
import io
import os
import random
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData

    set_gene = set()

    print 'reading from... {}'.format(path_inputData+'/pathway_cancer.xml')
    k = 0
    entity_type = None
    gene_name = None
    f = open(path_inputData+'/pathway_cancer.xml', 'r')
    for line in f:
        k += 1
        if k % 5 == 1:
            line = line.strip().split(' ')[-1]
            entity_type = line.split('\"')[1]
        if k % 5 == 3:
            if entity_type != 'gene': continue
            line = line.strip().split('\"')[1]
            line = line.split(',')[0]
            gene_name = line.split('.')[0]
            gene_name = gene_name.lower()
            set_gene.add(gene_name)

    f = open(path_outputData+'/pathway_gene.txt', 'w')
    for line in set_gene:
        print >> f, line
    f.close()
#EOF.
