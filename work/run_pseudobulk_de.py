from pathlib import Path
import os,sys,json,gc,inspect
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np,pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
ROOT=Path(__file__).resolve().parent
def process(acc,configs):
    folder=ROOT/'processed'/acc;out=ROOT/'results'/acc;out.mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(folder/'sample_cell_counts.tsv',sep='\t')
    for config in configs:
        for typ in ['Astrocyte','Endothelial']:
            dest=out/f'{config}__{typ}__de.tsv.gz'
            if dest.exists() and not __import__('os').environ.get('ISCHEMIA_FORCE_RECOMPUTE'):print('EXISTS',dest.name,flush=True);continue
            sub=manifest[(manifest.config==config)&(manifest.cell_type==typ)&(manifest.n_cells>=30)].copy().set_index('sample')
            if sub.condition.value_counts().min()<2:print('INSUFFICIENT',config,typ,flush=True);continue
            cols=[]
            for sample in sub.index:
                x=pd.read_csv(folder/(sample+'_pseudobulk.tsv.gz'),sep='\t',index_col=0)[config+'__'+typ].rename(sample);cols.append(x)
            counts=pd.concat(cols,axis=1).fillna(0).astype(np.int64).T
            minimum=int(sub.condition.value_counts().min());keep=(counts>=10).sum(axis=0)>=minimum
            counts=counts.loc[:,keep]
            metadata=sub.loc[counts.index,['condition']].copy();metadata['condition']=pd.Categorical(metadata.condition,categories=['Sham','MCAO'])
            print('DE',acc,config,typ,'samples',len(counts),'genes',counts.shape[1],flush=True)
            dds=DeseqDataSet(counts=counts,metadata=metadata,design='~condition',refit_cooks=True,n_cpus=1,quiet=True)
            dds.deseq2();stats=DeseqStats(dds,contrast=['condition','MCAO','Sham'],n_cpus=1,quiet=True)
            stats.summary();result=stats.results_df
            result.to_csv(dest,sep='\t',index_label='gene')
            norm=pd.DataFrame(dds.layers['normed_counts'],index=counts.index,columns=counts.columns)
            norm.T.to_csv(out/f'{config}__{typ}__normalized.tsv.gz',sep='\t',index_label='gene')
            diag={'dataset':acc,'config':config,'cell_type':typ,'n_samples':len(counts),'samples_by_condition':metadata.condition.value_counts().to_dict(),'tested_genes':len(result),'fdr05':int((result.padj<.05).sum()),'fdr05_abslog2fc05':int(((result.padj<.05)&(result.log2FoldChange.abs()>=.5)).sum()),'min_cells':int(sub.n_cells.min()),'n_cpus':1,'software':'PyDESeq2','replication_level':'exploratory: validation MCAO n=2' if minimum<3 else 'sample-level'}
            (out/f'{config}__{typ}__diagnostics.json').write_text(json.dumps(diag,indent=2),encoding='utf-8');print(diag,flush=True)
            del dds,stats,counts,norm;gc.collect()
if __name__=='__main__':process(sys.argv[1] if len(sys.argv)>1 else 'GSE174574',sys.argv[2:] or ['primary'])
