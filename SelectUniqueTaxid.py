import argparse
    
def read_blast_table(blastx_file):
    '''Read blastx table'''
    totalList=[]

    for line in blastx_file:
        
        line=line.rstrip()
        if line.startswith('ID'):
            continue
        else:
            parts=line.split('\t')

            try:
                taxid=parts[int(opts.position)].split(';')
            except:
                taxid=''
            
            for t in taxid:
                if t != '.':
                    totalList.append(t)

    for t in sorted(set(totalList)):
        print(t)


def main():
    parser = argparse.ArgumentParser('Select unique taxid')
    parser.add_argument('-i_Blastx', '--fileBlast', help="blastx output table")
    parser.add_argument('-position', '--position', help="column position of taxid")

    global opts
    opts = parser.parse_args()

    in_fileBlastx = open(opts.fileBlast)
    read_blast_table(in_fileBlastx)

main()