"""Post-QC amendment: remove cells with a strong mixed myeloid signature.
This deliberately conservative sensitivity is not an ambient-RNA correction.
"""
from pathlib import Path
import sys,gc
import numpy as np,pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parent
def process(acc):
 p=ROOT/'processed'/acc;records=[]
 for f in sorted(p.glob('*_target_counts.npz')):
  s=f.name.split('_')[0];x=sparse.load_npz(f);q=pd.read_csv(p/f'{s}_target_cells.tsv.gz',sep='\t');g=pd.Index(pd.read_csv(p/f'{s}_target_genes.tsv.gz',sep='\t').symbol)
  panel=['Ptprc','Tyrobp','Lyz2','C1qa','C1qb','C1qc','Csf1r','Aif1'];y=x[:,g.get_indexer(panel)];cs=np.asarray(y.sum(axis=1)).ravel();nd=np.asarray((y>0).sum(axis=1)).ravel();cd45=x[:,g.get_loc('Ptprc')].toarray().ravel()
  flag=((cs/np.maximum(q.n_umis.to_numpy(),1)>.02)&(nd>=3))|(cd45>=3)
  q['mixed_myeloid_flag']=flag;q['myeloid_umi_fraction']=cs/q.n_umis
  q.to_csv(p/f'{s}_target_cells.tsv.gz',sep='\t',index=False)
  condition='Sham' if ('sham' in q.prefix.iloc[0].lower() or '_WTC' in q.prefix.iloc[0]) else 'MCAO'
  pb=pd.read_csv(p/f'{s}_pseudobulk.tsv.gz',sep='\t');fr=pd.read_csv(p/f'{s}_fractions.tsv.gz',sep='\t');me=pd.read_csv(p/f'{s}_mean_logcp10k.tsv.gz',sep='\t')
  assert (pb.gene==g).all() and (fr.gene==g).all() and (me.gene==g).all()
  for typ in ['Astrocyte','Endothelial','Pericyte']:
   idx=np.flatnonzero(q.cell_type.eq(typ).to_numpy()&~flag);y=x[idx];n=len(idx);key='low_myeloid__'+typ
   pb[key]=np.asarray(y.sum(axis=0)).ravel();fr[key]=np.asarray((y>0).sum(axis=0)).ravel()/n
   y=y.astype(float).multiply((1e4/q.iloc[idx].n_umis.to_numpy())[:,None]).tocsr();y.data=np.log1p(y.data);me[key]=np.asarray(y.mean(axis=0)).ravel()
   records.append(dict(dataset=acc,sample=s,condition=condition,config='low_myeloid',cell_type=typ,n_cells=n,total_target_counts=int(pb[key].sum()),removed_cells=int((q.cell_type.eq(typ)&flag).sum())))
  pb.to_csv(p/f'{s}_pseudobulk.tsv.gz',sep='\t',index=False);fr.to_csv(p/f'{s}_fractions.tsv.gz',sep='\t',index=False);me.to_csv(p/f'{s}_mean_logcp10k.tsv.gz',sep='\t',index=False)
  del x,y;gc.collect()
 manifest=pd.read_csv(p/'sample_cell_counts.tsv',sep='\t');manifest=manifest[manifest.config!='low_myeloid'];r=pd.DataFrame(records);pd.concat([manifest,r],ignore_index=True).to_csv(p/'sample_cell_counts.tsv',sep='\t',index=False)
 print(r.to_string(index=False))
if __name__=='__main__':process(sys.argv[1])
