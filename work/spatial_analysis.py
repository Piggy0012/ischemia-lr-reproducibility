"""Descriptive expression maps using sample-specific deposited alignment coordinates.
No image-defined anatomy is inferred. Spots are not treated as biological replicates.
"""
from pathlib import Path
import gzip,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from stream_annotate import chunks
ROOT=Path(__file__).resolve().parent;RAW=ROOT/'data/GSE233814/raw';OUT=ROOT/'results/spatial';OUT.mkdir(exist_ok=True)
geometry=pd.read_csv(ROOT/'spatial_geometry/visium_row_col.tsv',sep='\t')
# SpatialFeatureExperiment's R lookup is 1-based; Loupe JSON is 0-based.
assert geometry.row.min()==1 and geometry.row.max()==78 and geometry.col.min()==1 and geometry.col.max()==128
geometry[['row','col']]=geometry[['row','col']]-1
panels=json.loads((ROOT/'gene_panels.json').read_text())
extra=['Aldh1l1','Slc1a2','Slc1a3','Glul','Pecam1','Cdh5','Kdr','Spp1','Itga5','Itgb1','Ptn','Ptprb','Ptprz1','Plat','Lrp1','Timp3','S1pr1','Timp1','Cd63','Tyrobp','C1qa','Lyz2']
genes=sorted(set(extra+sum([v['genes'] for v in panels.values()],[])));rows=[];stats=[]
times={'GSM7437221':'Control','GSM7437222':'Day 1','GSM7437223':'Day 3','GSM7437224':'Day 7a','GSM7437225':'Day 7b'}
for f in sorted(RAW.glob('*_matrix.mtx.gz')):
 prefix=f.name.removesuffix('_matrix.mtx.gz');s=prefix.split('_')[0]
 fs=pd.read_csv(RAW/f'{prefix}_features.tsv.gz',sep='\t',header=None);bc=pd.read_csv(RAW/f'{prefix}_barcodes.tsv.gz',sep='\t',header=None)[0]
 mapping={g:i for i,g in enumerate(genes)};gi=np.array([mapping.get(g,-1) for g in fs[1]])
 it=chunks(f);ng,nc,nnz=next(it);assert nc==len(bc) and ng==len(fs)
 x=np.zeros((nc,len(genes)),np.int32);total=np.zeros(nc,np.float64);ngenes=np.zeros(nc,np.int32);seen=0
 for a in it:
  g,c,v=a.T;total+=np.bincount(c,weights=v,minlength=nc);ngenes+=np.bincount(c,minlength=nc).astype(np.int32);yes=gi[g]>=0;np.add.at(x,(c[yes],gi[g[yes]]),v[yes]);seen+=len(a)
 assert seen==nnz
 align=json.load(gzip.open(next(RAW.glob(s+'*.json.gz')),'rt'));points=pd.DataFrame(align['oligo'])
 assert len(geometry.merge(points,on=['row','col'],validate='one_to_one'))==4992
 q=pd.DataFrame({'barcode':bc,'barcode_core':bc.str.replace(r'-\d+$','',regex=True),'n_umis':total.astype(int),'n_genes':ngenes})
 q=q.merge(geometry.rename(columns={'barcode':'barcode_core'}),on='barcode_core',how='left',validate='one_to_one',sort=False).merge(points[['row','col','imageX','imageY']],on=['row','col'],how='left',validate='one_to_one',sort=False)
 assert q.imageX.notna().all() and q.imageY.notna().all() and len(q)==nc
 q['sample']=s;q['time']=times[s];q['qc_pass']=(q.n_genes>=200)&(q.n_umis>=500)
 norm=np.log1p(x/np.maximum(total[:,None],1)*1e4)
 for i,g in enumerate(genes):q[g]=norm[:,i]
 q['Astro_marker_mean']=q[['Aldh1l1','Slc1a2','Slc1a3']].mean(axis=1);q['EC_marker_mean']=q[['Pecam1','Cdh5','Kdr']].mean(axis=1)
 for name,panel in panels.items():q[name]=q[panel['genes']].mean(axis=1)
 q['Spp1_integrin_coexpression']=np.sqrt(q.Spp1*q[['Itga5','Itgb1']].min(axis=1));q['Ptn_Ptprb_coexpression']=np.sqrt(q.Ptn*q.Ptprb)
 q.to_csv(OUT/f'{s}_spots.tsv.gz',sep='\t',index=False)
 good=q[q.qc_pass].copy();coords=good[['imageX','imageY']].to_numpy();nn=cKDTree(coords).query(coords,k=7)[1][:,1:]
 for l,r in [('Spp1','Itga5'),('Ptn','Ptprz1'),('Plat','Lrp1'),('Timp3','Kdr')]:
  stats.append(dict(sample=s,time=times[s],ligand=l,receptor=r,n_spots=len(good),spearman_same_spot=float(spearmanr(good[l],good[r]).statistic),spearman_neighbor=float(spearmanr(good[l].to_numpy(),good[r].to_numpy()[nn].mean(axis=1)).statistic),ligand_fraction=float((good[l]>0).mean()),receptor_fraction=float((good[r]>0).mean())))
 rows.append(dict(sample=s,time=times[s],submitted_spots=nc,qc_spots=len(good),matched_coordinates=int(q.imageX.notna().sum()),min_umis=int(total.min()),median_umis=float(np.median(total)),coordinate_source='sample-specific alignment JSON joined by verified Visium barcode-to-row/col lookup'))
 print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(OUT/'sample_qc.tsv',sep='\t',index=False);pd.DataFrame(stats).to_csv(OUT/'descriptive_coexpression.tsv',sep='\t',index=False)
