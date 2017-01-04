from collections import defaultdict as dd
import io
import os
import random
import argparse

'''
python analysis2.py --inputData /usr1/public/yifeng/Github/inputData --outputData /usr1/public/yifeng/Github/inputData --pathwayName cancer
'''
def extract_xml(name_xml, name_pathway, path_inputData, path_outputData, all_gene_set):
    
    k = 0
    entity_type = None
    gene_name = None
    set_gene = set()
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
    print '{}\t{}\t{}/{}'.format(name_xml, name_pathway, len(all_gene_set.intersection(set_gene)), len(set_gene))
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

    all_gene_set = [line.strip() for line in open(path_outputData+'/DEGl.txt')]
    all_gene_set = set(all_gene_set)

    # pathway2xml = {'pi3k':'hsa04151','cancer':'hsa05200'}
    pathway2xml = {'coad':'hsa05210',
    'gbm':'hsa05214',
    'blca':'hsa05219',
    'prad': 'hsa05215',
    'ucec': 'hsa05213',
    'brca': 'hsa05224',
    'pi3k':'hsa04151',
    'cancer':'hsa05200',
    'p53': 'hsa04115',
    'notch': 'hsa04330'}
    print 'KEGGpathway\tcancer\tintersect/pathway'
    for name_pathway in pathwayName:
        name_xml = pathway2xml[name_pathway]
        extract_xml(name_xml, name_pathway, path_inputData, path_outputData, all_gene_set)



#EOF.
