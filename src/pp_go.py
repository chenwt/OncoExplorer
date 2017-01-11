# process the go term files.

if __name__ == '__main__':

    # Process goAnnotation.
    goAnn = set()
    path_goAnn = '/usr1/public/yifeng/Github/inputData/goa_human.gaf'
    f = open(path_goAnn, 'r')
    for line in f:
        if line[0] == '!': continue
        line = line.split('\t')
        gene, goId = line[2], line[4]
        goAnn.add((gene.lower(), goId.lower()))
        #print gene, goId
        #print line

    f.close()

    path_goAnn = '/usr1/public/yifeng/Github/PropprData/goAnn_full.cfacts'
    f = open(path_goAnn, 'w')
    for line in goAnn:
        print >> f, 'goAnn\t'+line[0]+'\tgo'+line[1].split(':')[1]
    f.close()




    # Getting the isA.cfacts file.
    path_inputData = '/usr1/public/yifeng/Github/inputData'
    path_outputData = '/usr1/public/yifeng/Github/PropprData'



    f = open(path_inputData+'/go-basic.obo', 'r')
    # Format: isA GoTermX GoTermY
    # Meaning: GoTermX is GoTermY
    fo = open(path_outputData+'/isA.cfacts', 'w')


    go_bp_set = set() # The set that belongs to biological process.

    flag_term = False
    #flag_biopro = true
    src_id = ''
    name_space = ''
    dst_id = list()
    for line in f:
        line = line.strip()
        if line == '[Term]': flag_term = True

        if flag_term == False: continue

        # Switch...case...
        if line.startswith('id: GO:'):
            src_id = line
            continue
        if line.startswith('namespace:'):
            name_space = line
            continue
        if line.startswith('is_a: GO:'):
            line = line.split(':')[2].split(' ! ')[0]
            dst_id.append('go'+line)
            continue
        if line == '':
            
            if name_space == 'namespace: biological_process':
            #print >> fo, '======================='
                src = 'go'+src_id.split(':')[2]

                go_bp_set.add(src)

                for dst in dst_id:
                    go_bp_set.add(dst)
                    print >> fo, 'isA\t'+src+'\t'+dst
            dst_id = list()
            flag_term = False
            continue
        #print line
    f.close()
    fo.close()


    path_goAnn = '/usr1/public/yifeng/Github/PropprData/goAnn.cfacts'
    f = open(path_goAnn, 'w')
    for line in goAnn:
        goterm = 'go'+line[1].split(':')[1]
        if goterm in go_bp_set:
            print >> f, 'goAnn\t'+line[0]+'\t'+goterm
    f.close()


    # Processing test.examples.
    cancer = 'ov'
    path_lut = '/usr1/public/yifeng/Github/outputData/lut_'+cancer+'.txt'
    f = open(path_lut, 'r')
    next(f)

    gene_set = set()
    for line in f:
        gene = line.strip().split('\t')[0]
        gene_set.add(gene)



    f.close()

    path_exam = '/usr1/public/yifeng/Github/PropprData/test_'+cancer+'.examples'
    f = open(path_exam, 'w')
    for line in gene_set:
        print >> f, 'isBP('+line+',Y)'

    f.close()

    base_set = set()
    path_base = '/usr1/public/yifeng/Github/PropprData/isBase.bak'
    f = open(path_base, 'r')
    for line in f:
        base = line.strip().split('\t')[1]
        base_set.add(base)
    f.close()

    path_base = '/usr1/public/yifeng/Github/PropprData/isBase.cfacts'
    f = open(path_base, 'w')
    for line in base_set:
        print >> f, 'isBase\t'+line
    f.close()

#EOF.
