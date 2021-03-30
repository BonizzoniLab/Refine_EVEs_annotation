import argparse
import re

class info():
    ref=''
    alt=''
    pass

class blastx_hit():
    q_seq_id=''
    q_start=''
    q_end=''
    s_id=''
    accession=''
    evalue=-1
    q_frame=-100
    perc_id=-1
    q_perc_cov=-1
    s_start=''
    s_end=''
    s_len=''
    taxid=''
    pass

class subject_summary():
    num_total_hit=0
    num_viral_hit=0
    hit_db={}
    pass

class Merged_region():
    scaffold=''
    merged_start=''
    merged_end=''
    topHit_name=''
    bestEval_name=''
    id=''
    pass

def read_VHname(VHfile):
    '''Read VH-classifier file virus.csv'''
    
    vtaxid=[]

    for line in VHfile:
        
        line=line.rstrip()
        parts=line.split(',')
        vtaxid.append(parts[0])
    
    return sorted(vtaxid)

def read_bedTopHits(in_fileTopHits):
    '''reads the top hits bed file of whitefield 2017'''
    db_th={}

    i=1
    for line in in_fileTopHits:
        
        line=line.rstrip()
        parts=line.split('\t')
        merged_reg=Merged_region()
        merged_reg.scaffold=parts[0]
        merged_reg.merged_start=int(parts[1])
        merged_reg.merged_end=int(parts[2])
        merged_reg.topHit_name=str(parts[0])+":"+str(parts[1])+"-"+str(parts[2])
        #merged_reg.topHit_name=parts[3]
        merged_reg.id='NIRVS_'+str(i)
        
        
        db_th[merged_reg.topHit_name]=merged_reg
        i=i+1
        
    return db_th

def read_blast_table(blastx_file):
    '''Read blastx table'''
    db={}
    i=1

    for line in blastx_file:
        
        line=line.rstrip()
        parts=line.split('\t')
        instance=blastx_hit()
        instance.q_seq_id=parts[0]
        instance.q_start=parts[1]
        instance.q_end=parts[2]
        instance.s_id=parts[3]
        instance.accession=parts[4]
        instance.evalue=float(parts[5])
        instance.q_frame=int(parts[6])
        instance.perc_id=float(parts[7])
        instance.q_perc_cov=float(parts[8])
        instance.s_start=parts[9]
        instance.s_end=parts[10]
        instance.s_len=parts[11]
        instance.taxid=parts[12]
        
        db[i]=instance
        i=i+1
    
    return db

def best_data_selection(db, db_topHits,  vtaxid_list):
    '''Selection of the best data between hits with the same query sequence id'''
    db_blastx={}
    complete_table=open(opts.output+'/CompleteTable.txt','w')

    for query in sorted(db_topHits.keys()): # in query_list:
        summary_ist=subject_summary()
        hit_db={}
        
        for v in db.itervalues():
            check_at_least_one_viral=0
            if v.q_seq_id == query:
                summary_ist.num_total_hit=summary_ist.num_total_hit+1
                if int(v.s_start) > int(v.s_end):
                    hit_db[str(v.accession)+'_'+str(v.s_id)+':'+str(v.s_end)+'-'+str(v.s_start)]=v
                else:
                    hit_db[str(v.accession)+'_'+str(v.s_id)+':'+str(v.s_start)+'-'+str(v.s_end)]=v
                
                # v_taxid_list=[]
                for taxid in v.taxid.split(';'):
                    if taxid in vtaxid_list:
                        check_at_least_one_viral=1
                #         if taxid not in v_taxid_list:
                #             v_taxid_list.append(taxid)
                
                # v_taxid_str=';'.join(v_taxid_list)
                # v.taxid=v_taxid_str

                if int(check_at_least_one_viral) == 1:
                    summary_ist.num_viral_hit=summary_ist.num_viral_hit+1

        summary_ist.hit_db=hit_db

        if int(summary_ist.num_viral_hit) > int(summary_ist.num_total_hit):
            print(str(query)+'\t'+str(summary_ist.num_total_hit)+'\t'+str(summary_ist.num_viral_hit))  

        db_blastx[query]=summary_ist
               
    complete_table.write('ID\tScaffold\tStart\tEnd\tLen\tq_seq_id\tNum_total_hit\tNum_viral_hit\tFraction_viral_hit\t' \
        'Bestviral_Protein\tVirus_name\tBestviralEval\tNum_subj_same_bestViralEval\t' \
        'Bestviral_frame\tBestviral_percid\tBestviral_qcov\tBestviral_start\tBestviral_end\tBestviral_slen\tBestviral_taxid\tBestviral_accession\t'
        'BestSubj\tBestEval\tNum_subj_same_bestEval\t' \
        'Sel_BestSubj\tSel_BestEval\tSel_Num_subj_same_bestEval\t' \
        'NIRVSmerged\n')

    for k,v in db_blastx.iteritems():
        q_covs=[]
        nonviral_evalues=[]
        nonviral_elem=[]
        selected_nonviral_evalues=[]
        selected_nonviral_elem=[]
        viral_evalues=[]
        viral_elem=[]

        min_viral_items=[]
        viral_items_strand='.'
        viral_items_identity='.'
        viral_items_qcov='.'
        viral_items_start='.'
        viral_items_end='.'
        viral_items_slen='.'
        viral_items_taxid='.'
        viral_items_accession='.'
        Bestviral_Protein=''
        Virus_name=''

        min_nonviral_items=[]
        min_nonviral_items_start=''
        min_nonviral_items_end=''

        selected_min_nonviral_items=[]
        selected_min_nonviral_items_start=''
        selected_min_nonviral_items_end=''

        min_nonviral_eval=-1
        selected_min_nonviral_eval=-1
        min_viral_eval=-1
        max_q_cov=-1

        Num_subj_same_bestViralEval=0
        Num_subj_same_bestNonViralEval=0
        selected_Num_subj_same_bestNonViralEval=0
        fragmented='no'

        for k_hit_db, v_hit_db in v.hit_db.iteritems():
            for taxid in v_hit_db.taxid.split(';'):
                if taxid in vtaxid_list:
                    q_covs.append(v_hit_db.q_perc_cov)

                    viral_evalues.append(v_hit_db.evalue)
                    viral_elem.append(k_hit_db)
                else:
                    nonviral_evalues.append(v_hit_db.evalue)
                    nonviral_elem.append(k_hit_db)
                    
                    if re.search(r'.*AGAP.*-PA.*', k_hit_db) is None and \
                        re.search(r'.*AAEL.*-P.*', k_hit_db) is None and \
                        re.search(r'.*ncharacterized.*', k_hit_db) is None and \
                        re.search(r'.*PREDICTED.*', k_hit_db) is None and \
                        re.search(r'.*hypothetical.*', k_hit_db) is None and \
                        re.search(r'.*CLUMA_CG.*', k_hit_db) is None and \
                        re.search(r'.*baculoviral.*', k_hit_db) is None and \
                        re.search(r'.*LOW QUALITY PROTEIN:.*', k_hit_db) is None and \
                        re.search(r'.*unnamed protein product.*', k_hit_db) is None:
                        selected_nonviral_evalues.append(v_hit_db.evalue)
                        selected_nonviral_elem.append(k_hit_db)

        if q_covs !=[]:
            max_q_cov=max(q_covs)

        if viral_evalues !=[]:
            indices=[]
            sel_items=[]
            pid=[]
            items=[]
            min_viral_eval=min(viral_evalues)
            indices = [i for i, x in enumerate(viral_evalues) if float(x) == float(min_viral_eval)]
            #min_viral_items=viral_elem[viral_evalues.index(min_viral_eval)]
            if len(indices) > int(1):
                for ind in indices:
                    items.append(viral_elem[ind])
                
                for i in items:
                    if float(v.hit_db[i].q_perc_cov) > float(80):
                        sel_items.append(i)
                        pid.append(v.hit_db[i].perc_id)
                
                if pid !=[]:
                    max_pid=max(pid)
                    max_pid_viral_item=sel_items[pid.index(max_pid)]
                    min_viral_items=max_pid_viral_item
                else:
                    min_viral_items=viral_elem[viral_evalues.index(min_viral_eval)]

            else:
                min_viral_items=viral_elem[viral_evalues.index(min_viral_eval)]

            Bestviral_Protein=min_viral_items.split('[')[0]
            Virus_name=min_viral_items.split('[')[1].split(']')[0]
            viral_items_strand=v.hit_db[min_viral_items].q_frame
            viral_items_identity=v.hit_db[min_viral_items].perc_id
            viral_items_qcov=v.hit_db[min_viral_items].q_perc_cov
            viral_items_start=v.hit_db[min_viral_items].s_start
            viral_items_end=v.hit_db[min_viral_items].s_end
            viral_items_slen=v.hit_db[min_viral_items].s_len
            viral_items_taxid=v.hit_db[min_viral_items].taxid
            viral_items_accession=v.hit_db[min_viral_items].accession
            Num_subj_same_bestViralEval=viral_evalues.count(min_viral_eval)
        else:
            Bestviral_Protein='NA'
            Virus_name='NA'

        if nonviral_evalues !=[]:
            min_nonviral_eval=min(nonviral_evalues)
            min_nonviral_items=nonviral_elem[nonviral_evalues.index(min_nonviral_eval)]
            Num_subj_same_bestNonViralEval=nonviral_evalues.count(min_nonviral_eval)
            
            if selected_nonviral_evalues !=[]:
                selected_min_nonviral_eval=min(selected_nonviral_evalues)
                selected_min_nonviral_items=selected_nonviral_elem[selected_nonviral_evalues.index(selected_min_nonviral_eval)]
                selected_Num_subj_same_bestNonViralEval=selected_nonviral_evalues.count(selected_min_nonviral_eval)
            else:
                selected_min_nonviral_items='NA'

        else:
            min_nonviral_items='NA'
            selected_min_nonviral_items='NA'
        
        if min_viral_items != 'NA' and float(max_q_cov) <= float(80):
            fragmented='yes'

        try:
            fraction_viral_hit=float(v.num_viral_hit)/float(v.num_total_hit)
        except:
            fraction_viral_hit=0

        complete_table.write(str(db_topHits[k].id)+'\t'+str(db_topHits[k].scaffold)+'\t'+str(db_topHits[k].merged_start)+'\t' \
            +str(db_topHits[k].merged_end)+'\t'+str(int(db_topHits[k].merged_end)-int(db_topHits[k].merged_start))+'\t' \
            +str(k)+'\t'+str(v.num_total_hit)+'\t'+str(v.num_viral_hit)+'\t'+str(fraction_viral_hit)+'\t' \
            +str(Bestviral_Protein)+'\t'+str(Virus_name)+'\t'+str(min_viral_eval)+'\t'+str(Num_subj_same_bestViralEval)+'\t' \
            +str(viral_items_strand)+'\t'+str(viral_items_identity)+'\t'+str(viral_items_qcov)+'\t' \
            +str(viral_items_start)+'\t'+str(viral_items_end)+'\t'+str(viral_items_slen)+'\t'+str(viral_items_taxid)+'\t'+str(viral_items_accession)+'\t' \
            +str(min_nonviral_items)+'\t'+str(min_nonviral_eval)+'\t'+str(Num_subj_same_bestNonViralEval)+'\t' \
            +str(selected_min_nonviral_items)+'\t'+str(selected_min_nonviral_eval)+'\t'+str(selected_Num_subj_same_bestNonViralEval)+'\t' \
            +str(fragmented)+'\n')
        
    complete_table.close()


def main():
    parser = argparse.ArgumentParser('Select viral results')
    parser.add_argument('-i_VHname', '--fileVH', help="tab delimited file")
    parser.add_argument('-i_Blastx', '--fileBlast', help="blastx output table")
    parser.add_argument('-i_TopHits', '--fileBedTopHits', help="bed file of the top hits results, Whitefield 2017")
    parser.add_argument('-output_path', '--output', help="path for the output files")
    
    global opts
    opts = parser.parse_args()

    in_fileVH = open(opts.fileVH)
    vtaxid_list=read_VHname(in_fileVH)

    in_fileTopHits = open(opts.fileBedTopHits)
    db_topHits=read_bedTopHits(in_fileTopHits)

    in_fileBlastx = open(opts.fileBlast)
    db=read_blast_table(in_fileBlastx)

    best_data_selection(db, db_topHits, vtaxid_list)
    

main()