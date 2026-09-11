from pathlib import Path
import json,sys,hashlib,importlib.metadata
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parent;checks=[]
def check(name,condition,details=None):
 assert bool(condition),name
 checks.append({'check':name,'status':'PASS','details':details})
for acc in ['GSE174574','GSE245386']:
 p=ROOT/'processed'/acc;m=pd.read_csv(p/'sample_cell_counts.tsv',sep='\t')
 check(acc+' unique sample/config/type',not m.duplicated(['sample','config','cell_type']).any())
 for s in m['sample'].unique():
  x=pd.read_csv(p/f'{s}_pseudobulk.tsv.gz',sep='\t',index_col=0);fr=pd.read_csv(p/f'{s}_fractions.tsv.gz',sep='\t',index_col=0)
  check(s+' nonnegative integer unique genes',x.index.is_unique and (x.to_numpy()>=0).all() and np.allclose(x,x.round()))
  check(s+' valid detection fractions',fr.index.equals(x.index) and (fr.to_numpy()>=0).all() and (fr.to_numpy()<=1).all())
  for r in m[m['sample']==s].itertuples(index=False):check(s+' count conservation '+r.config+' '+r.cell_type,int(x[r.config+'__'+r.cell_type].sum())==int(r.total_target_counts))
 q=pd.read_csv(p/'qc_summary.tsv',sep='\t');check(acc+' cell accounting',((q.n_qc<=q.n_input)&(q.n_assigned<=q.n_qc)).all())
 c=pd.read_csv(ROOT/f'results/{acc}/primary__communication.tsv.gz',sep='\t');sc=pd.read_csv(ROOT/f'results/{acc}/primary__communication_sample_scores.tsv.gz',sep='\t')
 check(acc+' complex score arithmetic',np.allclose(sc.score,np.sqrt(sc.ligand_logexpr*sc.receptor_logexpr),rtol=1e-10))
 check(acc+' exact p range',(c.exact_permutation_p>=.1-1e-12).all() and (c.exact_permutation_p<=1).all())
 check(acc+' nested detection gates',((~c.eligible_20)|c.eligible_10).all() and ((~c.eligible_10)|c.eligible_05).all())
 for config in ['primary','strict_identity','mt10','low_myeloid']:
  for typ in ['Astrocyte','Endothelial']:
   d=pd.read_csv(ROOT/f'results/{acc}/{config}__{typ}__de.tsv.gz',sep='\t');diag=json.loads((ROOT/f'results/{acc}/{config}__{typ}__diagnostics.json').read_text())
   check(acc+config+typ+' DE count',len(d)==diag['tested_genes'] and int((d.padj<.05).sum())==diag['fdr05'])
s=json.loads((ROOT/'results/cross_cohort/summary.json').read_text());c=pd.read_csv(ROOT/'results/cross_cohort/communication_replication.tsv.gz',sep='\t');common=c[c.both_eligible10]
check('cross-cohort 323 / 258',len(common)==323 and common.same_direction.sum()==258)
selected=common[common.discovery_selected];check('discovery selected 182 / 154',len(selected)==182 and selected.same_direction.sum()==154)
check('correlation',np.isclose(spearmanr(common.score_difference_discovery,common.score_difference_validation).statistic,s['communication']['rho_all_common']))
sp=pd.read_csv(ROOT/'results/spatial/sample_qc.tsv',sep='\t');check('spatial spots and geometry',sp.qc_spots.sum()==11969 and sp.submitted_spots.sum()==12099 and (sp.matched_coordinates==sp.submitted_spots).all())
for f in (ROOT/'results/spatial').glob('*_spots.tsv.gz'):
 d=pd.read_csv(f,sep='\t');check(f.name+' coordinate integrity',d.barcode.is_unique and d.imageX.notna().all() and d.imageY.notna().all() and not d.duplicated(['row','col']).any())
refs=json.loads((ROOT.parent/'outputs/references.json').read_text(encoding='utf-8'));check('22 cited references verified',len(refs)==22 and all(r['metadata_status']=='crossref_verified' for r in refs))
pkgs=['numpy','pandas','scipy','statsmodels','pydeseq2','anndata','h5py','matplotlib','seaborn','scikit-learn','requests','rdata','python-docx']
versions={p:importlib.metadata.version(p) for p in pkgs};versions['python']=sys.version
(ROOT/'software_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
(ROOT/'requirements_analysis.txt').write_text('\n'.join(f'{p}=={versions[p]}' for p in pkgs)+'\n',encoding='utf-8')
(ROOT/'results/validation_checks.json').write_text(json.dumps({'checks':checks,'n_passed':len(checks),'scope':'arithmetic, mappings, result consistency; not independent biological replication'},indent=2),encoding='utf-8')
print('PASS',len(checks),'checks; versions',versions)
