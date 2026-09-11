from pathlib import Path
import sys,gzip,json,gc
import numpy as np,pandas as pd
from scipy import sparse
from stream_annotate import chunks
ROOT=Path(__file__).resolve().parent
def process(acc):
    raw=ROOT/'data'/acc/'raw';out=ROOT/'processed'/acc
    for matrix in sorted(raw.glob('*matrix.mtx.gz')):
        prefix=matrix.name.removesuffix('_matrix.mtx.gz');sample=prefix.split('_')[0]
        dest=out/(sample+'_pseudobulk.tsv.gz')
        if dest.exists(): print('EXISTS',sample,flush=True);continue
        q=pd.read_csv(out/(sample+'_cells.tsv.gz'),sep='\t');features=pd.read_csv(out/(sample+'_features.tsv.gz'),sep='\t')
        wanted=q.cell_type.isin(['Astrocyte','Endothelial','Pericyte']).to_numpy();selected=np.flatnonzero(wanted)
        mapping=np.full(len(q),-1,np.int32);mapping[selected]=np.arange(len(selected))
        nnz=int(q.loc[wanted,'n_genes'].sum());ri=np.empty(nnz,np.int32);ci=np.empty(nnz,np.int32);val=np.empty(nnz,np.int32);pos=0
        print('AGGREGATE',sample,'cells',len(selected),'nnz',nnz,flush=True)
        iterator=chunks(matrix);ng,nc,nnz_all=next(iterator)
        for a in iterator:
            yes=wanted[a[:,1]];a=a[yes];n=len(a)
            ri[pos:pos+n]=mapping[a[:,1]];ci[pos:pos+n]=a[:,0];val[pos:pos+n]=a[:,2];pos+=n
        assert pos==nnz
        x=sparse.coo_matrix((val,(ri,ci)),shape=(len(selected),ng)).tocsr();del ri,ci,val
        # Collapse duplicate symbols before detecting genes or testing expression.
        symbols=features.symbol.astype(str);unique=pd.Index(symbols.unique());ids=unique.get_indexer(symbols)
        if len(unique)<ng:
            collapse=sparse.csr_matrix((np.ones(ng,np.int32),(np.arange(ng),ids)),shape=(ng,len(unique)))
            x=(x@collapse).tocsr()
        q=q.iloc[selected].copy().reset_index(drop=True)
        q.to_csv(out/(sample+'_target_cells.tsv.gz'),sep='\t',index=False)
        pd.DataFrame({'symbol':unique}).to_csv(out/(sample+'_target_genes.tsv.gz'),sep='\t',index=False)
        sparse.save_npz(out/(sample+'_target_counts.npz'),x)
        group='Sham' if ('sham' in prefix.lower() or '_WTC' in prefix) else 'MCAO'
        configs={'primary':np.ones(len(q),bool),'strict_identity':q.strict_identity.to_numpy(),'mt10':(q.pct_mt<=10).to_numpy()}
        records=[];pseudobulk={'gene':unique};fractions={'gene':unique};means={'gene':unique}
        for config,mask in configs.items():
            for typ in ['Astrocyte','Endothelial','Pericyte']:
                idx=np.flatnonzero(mask&(q.cell_type==typ).to_numpy());y=x[idx];n=len(idx);key=config+'__'+typ
                pseudobulk[key]=np.asarray(y.sum(axis=0)).ravel()
                if n:
                    fractions[key]=np.asarray((y>0).sum(axis=0)).ravel()/n
                    # Mean of per-cell log1p(CP10k), matching ligand/receptor detection gate.
                    y=y.astype(np.float64).multiply((1e4/np.maximum(q.iloc[idx].n_umis.to_numpy(),1))[:,None]).tocsr();y.data=np.log1p(y.data)
                    means[key]=np.asarray(y.mean(axis=0)).ravel()
                else:fractions[key]=np.zeros(len(unique));means[key]=np.zeros(len(unique))
                records.append({'dataset':acc,'sample':sample,'condition':group,'config':config,'cell_type':typ,'n_cells':n,'total_target_counts':int(np.sum(pseudobulk[key]))})
        pd.DataFrame(pseudobulk).to_csv(dest,sep='\t',index=False)
        pd.DataFrame(fractions).to_csv(out/(sample+'_fractions.tsv.gz'),sep='\t',index=False)
        pd.DataFrame(means).to_csv(out/(sample+'_mean_logcp10k.tsv.gz'),sep='\t',index=False)
        pd.DataFrame(records).to_csv(out/(sample+'_sample_counts.tsv'),sep='\t',index=False)
        print(records[:3],flush=True);del x,y;gc.collect()
    allrec=pd.concat([pd.read_csv(f,sep='\t') for f in out.glob('*_sample_counts.tsv')],ignore_index=True)
    allrec.to_csv(out/'sample_cell_counts.tsv',sep='\t',index=False)
if __name__=='__main__':process(sys.argv[1] if len(sys.argv)>1 else 'GSE174574')
