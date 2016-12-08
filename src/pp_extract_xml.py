from collections import defaultdict as dd
import io
import os
import random
import argparse

def extract_xml(name_xml, name_pathway, path_inputData, path_outputData):
    print 'reading from... {}'.format(path_inputData+'/'+name_xml+'.xml')
    k = 0
    entity_type = None
    gene_name = None
    f = open(path_inputData+'/'+name_xml+'.xml', 'r')
    for _ in range(7): next(f)
    for line in f:
        k += 1
        if k % 5 == 1:
            line = line.strip().split(' ')[-1]
            entity_type = line.split('\"')[1]
            if entity_type == 'group': break
        if k % 5 == 3:
            if entity_type != 'gene': continue
            line = line.strip().split('\"')[1]
            line = line.split(',')[0]
            gene_name = line.split('.')[0]
            gene_name = gene_name.lower()
            set_gene.add(gene_name)

    f = open(path_outputData+'/'+name_pathway+'.txt', 'w')
    for line in set_gene:
        print >> f, line
    f.close()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    parser.add_argument('--pathwayName', help = 'name of the pathway', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData
    pathwayName = args.pathwayName
    pathwayName = pathwayName.split(',')
    set_gene = set()

    pathway2xml = {'pi3k':'hsa04151','cancer':'hsa05200'}
    for name_pathway in pathwayName:
        name_xml = pathway2xml[name_pathway]
        extract_xml(name_xml, name_pathway, path_inputData, path_outputData)



#EOF.
