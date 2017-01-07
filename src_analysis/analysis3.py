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
    # Format: isA GoTermX GoTermY
    # Meaning: GoTermX is GoTermY
    fo = open(path_outputData+'/isA.cfacts', 'w')

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
            # if dst_id != '': dst_id += '\n' 
            # dst_id += line
            continue
        if line == '':
            
            if name_space == 'namespace: biological_process':
            #print >> fo, '======================='
                src = 'go'+src_id.split(':')[2]

                for dst in dst_id:
                    print >> fo, 'isA\t'+src+'\t'+dst
                # print >> fo, src
                # #print >> fo, name_space
                # print >> fo, ' '.join(dst_id)
                #for dst in dst_id:
                #    print >> fo, dst
                # if dst_id == '':
                #     dst_id = 'terminal'
                # else:
                #     dst_id = dst_id.split(':')[2]
                # print >> fo, dst_id
            dst_id = list()
            flag_term = False
            continue
        #print line
    f.close()
    fo.close()
#EOF.
