import argparse

if __name__ == '__main__':
    ''' 
    python analysis3.py --inputData /usr1/public/yifeng/Github/inputData --outputData /usr1/public/yifeng/Github/outputData
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputData', help = 'directory of input data', type = str)
    parser.add_argument('--outputData', help = 'directory of output data', type = str)
    args = parser.parse_args()

    path_inputData = args.inputData 
    path_outputData = args.outputData

    f = open(path_inputData+'/go-basic.obo', 'r')
    fo = open(path_outputData+'/process.txt', 'w')
    flag_term = False
    #flag_biopro = true
    src_id = ''
    name_space = ''
    dst_id = ''
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
            if dst_id != '': dst_id += '\n' 
            dst_id += line
            continue
        if line == '':
            
            if name_space == 'namespace: biological_process':
            #print >> fo, '======================='
                print >> fo, src_id
                print >> fo, name_space
                print >> fo, dst_id
            dst_id = ''
            flag_term = False
            continue
        #print line
    f.close()
    fo.close()
#EOF.
