from pathlib import Path
import json,gc
import numpy as np,pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parent
rows=[]
genes=['Aldh1l1','Slc1a2','Slc1a3','Glul','Sox9','Pecam1','Cdh5','Kdr','Ptprc','Tyrobp','Lyz2','C1qa','P2ry12','Tmem119','Spp1','Cd14','Timp1','Cd63','Ptn','Ptprb','S1pr1']
for acc in ['GSE174574','GSE245386']:
    p=ROOT/'processed'/acc
    for f in sorted(p.glob('*_target_counts.npz')):
        s=f.name.split('_')[0];x=sparse.load_npz(f)
        q=pd.read_csv(p/f'{s}_target_cells.tsv.gz',sep='\t');g=pd.Index(pd.read_csv(p/f'{s}_target_genes.tsv.gz',sep='\t').symbol)
        for typ in ['Astrocyte','Endothelial']:
            for config in ['primary','strict_identity']:
                mask=q.cell_type.eq(typ).to_numpy(copy=True)
                if config=='strict_identity':mask &=q.strict_identity.to_numpy()
                y=x[mask];num=int(mask.sum())
                for gene in genes:
                    if gene not in g:continue
                    v=y[:,g.get_loc(gene)].toarray().ravel()
                    rows.append(dict(dataset=acc,sample=s,cell_type=typ,config=config,gene=gene,n_cells=num,fraction=float((v>0).mean()),mean_umi=float(v.mean())))
        del x,y;gc.collect()
pd.DataFrame(rows).to_csv(ROOT/'results/identity_marker_audit.tsv',sep='\t',index=False)
d=pd.DataFrame(rows);d=d[d.config.eq('primary')]
print(d.groupby(['dataset','cell_type','gene']).fraction.mean().unstack('gene').round(3).to_string())
