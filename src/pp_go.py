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

    path_goAnn = '/usr1/public/yifeng/Github/inputData/goAnn.cfacts'
    f = open(path_goAnn, 'w')
    for line in goAnn:
        print >> f, line[0]+'\t'+line[1]
    f.close()

#EOF.
