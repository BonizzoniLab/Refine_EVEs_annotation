import argparse

class info():
    line=''
    taxid=''
    pass

class tax_info():
    OrderFamily=''
    OrderFamily_taxid=''
    pass

def read_fileCT(in_fileCT):
    '''reads output file of the selection of the best viral and not viral subject for each tophits of Whitefield et al 2017'''
    db_CT={}
    taxid=''

    for line in in_fileCT:
        
        line=line.rstrip()
        parts=line.split('\t')
        instance_info=info()
        instance_info.taxid=parts[19]
        instance_info.line=line
        db_CT[parts[0]]=instance_info
    
    return db_CT

def read_fileTK(in_fileTaxonkit):
    '''read output file of taxonkit reformat function in which viral order and family are selected for each viral taxid selected'''
    db_TK={}

    for line in in_fileTaxonkit:
        line=line.rstrip()
        parts=line.split('\t')
        inst_tk=tax_info()
        tk_taxid=parts[0]
        inst_tk.OrderFamily=parts[3]
        inst_tk.OrderFamily_taxid=parts[4]
        db_TK[tk_taxid]=inst_tk

    return db_TK

def AssignClass(db_CT, db_TK):
    '''Assign viral Order/family classification implemented with taxonkit according to NCBI taxonomy to the output classification'''
    out_file=open(str(opts.fileCT).split('.')[0]+'_classified.txt', 'w')

    for k in sorted(db_CT.keys()):
        v=db_CT[k]
        if k == 'ID':
            new_line=v.line+'\tViral Order\tViral Family\tViral Order taxid\tViral Family taxid\n'
            out_file.write(new_line)
        else:
            tax=v.taxid.split(';')[0]
            family='.'
            order='.'

            if tax in db_TK.keys():
                order=str(db_TK[tax].OrderFamily.split(';')[0])
                family=str(db_TK[tax].OrderFamily.split(';')[1])
                order_taxid=str(db_TK[tax].OrderFamily_taxid.split(';')[0])
                family_taxid=str(db_TK[tax].OrderFamily_taxid.split(';')[1])
                
                if family=='':
                    family='Undetermined'
                    family_taxid='Undetermined'
                if order=='':
                    order='Undetermined'
                    order_taxid='Undetermined'

                new_line=v.line+'\t'+order+'\t'+family+'\t'+order_taxid+'\t'+family_taxid+'\n'
                out_file.write(new_line)
            else:
                new_line=v.line+'\tNA\tNA\tNA\tNA\n'
                out_file.write(new_line)

    out_file.close()

def main():
    parser = argparse.ArgumentParser('Assign viral Order/family classification implemented with taxonkit according to NCBI taxonomy to the output classification')
    parser.add_argument('-CompleteTable', '--fileCT', help="CompleteTable_*.txt")
    parser.add_argument('-Taxonkit_out', '--fileTaxonkit', help="Taxonkit_reformat_Order-Family_Viral_taxid_*.txt")
    
    global opts
    opts = parser.parse_args()

    in_fileCT = open(opts.fileCT)
    db_CT=read_fileCT(in_fileCT)

    in_fileTaxonkit = open(opts.fileTaxonkit)
    db_TK=read_fileTK(in_fileTaxonkit)

    AssignClass(db_CT, db_TK)
    
main()