"""Aggregate original singlets and independently supported target cells."""
from pathlib import Path
import sys,json,gc,shutil,os
import numpy as np,pandas as pd
from revision_data import load_sample,samples
ROOT=Path(__file__).resolve().parent
ANALYSIS=ROOT/'revision_analysis'
CONFIGS=['singlet','reference_singlet','reference_only']

def process(acc,sample):
    out=ANALYSIS/'processed'/acc;out.mkdir(parents=True,exist_ok=True)
    dest=out/f'{sample}_pseudobulk.tsv.gz'
    if dest.exists() and os.environ.get('ISCHEMIA_REVISION_FORCE')!='1':print('AGGREGATE EXISTS',sample,flush=True);return
    x,q,genes=load_sample(acc,sample)
    r=pd.read_csv(ROOT/f'revision_results/identity/{acc}/{sample}_identity.tsv.gz',sep='\t')
    assert np.array_equal(q.barcode,r.barcode) and 'predicted_doublet' in r
    singlet=~r.predicted_doublet.to_numpy(bool);confident=r.whole_brain_probability.to_numpy()>=.5
    condition='Sham' if ('sham' in q.prefix.iloc[0].lower() or '_WTC' in q.prefix.iloc[0]) else 'MCAO'
    counts={'gene':genes};fraction={'gene':genes};mean={'gene':genes};records=[]
    for config in CONFIGS:
        for typ in ['Astrocyte','Endothelial','Pericyte']:
            original=(q.cell_type==typ).to_numpy();reference=(r.whole_brain_broad==typ).to_numpy()&confident
            mask=singlet&(original if config=='singlet' else (original&reference if config=='reference_singlet' else reference))
            idx=np.flatnonzero(mask);y=x[idx];n=len(idx);key=config+'__'+typ
            counts[key]=np.asarray(y.sum(axis=0)).ravel()
            if n:
                fraction[key]=np.asarray((y>0).sum(axis=0)).ravel()/n
                y=y.astype(np.float64).multiply((1e4/q.n_umis.to_numpy()[idx])[:,None]).tocsr();y.data=np.log1p(y.data)
                mean[key]=np.asarray(y.mean(axis=0)).ravel()
            else:fraction[key]=np.zeros(len(genes));mean[key]=np.zeros(len(genes))
            records.append({'dataset':acc,'sample':sample,'condition':condition,'config':config,'cell_type':typ,'n_cells':n,'total_target_counts':int(np.sum(counts[key]))})
    for suffix,values in [('pseudobulk',counts),('fractions',fraction),('mean_logcp10k',mean)]:
        pd.DataFrame(values).to_csv(out/f'{sample}_{suffix}.tsv.gz',sep='\t',index=False)
    pd.DataFrame(records).to_csv(out/f'{sample}_sample_counts.tsv',sep='\t',index=False)
    print('AGGREGATE DONE',sample,[(z['config'],z['cell_type'],z['n_cells']) for z in records],flush=True)
    del x,y;gc.collect()

if __name__=='__main__':
    for acc in sys.argv[1:] or ['GSE174574','GSE245386']:
        for sample in samples(acc):process(acc,sample)
        out=ANALYSIS/'processed'/acc
        pd.concat([pd.read_csv(f,sep='\t') for f in sorted(out.glob('*_sample_counts.tsv'))]).to_csv(out/'sample_cell_counts.tsv',sep='\t',index=False)
    (ANALYSIS/'literature').mkdir(parents=True,exist_ok=True)
    shutil.copy2(ROOT/'literature/mouseconsensus.csv',ANALYSIS/'literature/mouseconsensus.csv')
