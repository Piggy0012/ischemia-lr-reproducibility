from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results/integrated';OUT.mkdir(exist_ok=True)
PANELS=json.loads((ROOT/'gene_panels.json').read_text())
PAIRS=[('Astrocyte','Endothelial','Spp1','Itga5_Itgb1'),('Astrocyte','Endothelial','Timp3','Kdr'),('Endothelial','Astrocyte','Ptn','Ptprz1'),('Endothelial','Astrocyte','Plat','Lrp1')]
(OUT/'illustrative_pairs.json').write_text(json.dumps(PAIRS,indent=2),encoding='utf-8')
bbb=[];scores=[];components=[];pairrows=[];associations=[];sens=[]
for acc in ['GSE174574','GSE245386']:
 folder=ROOT/'processed'/acc;results=ROOT/'results'/acc
 manifest=pd.read_csv(folder/'sample_cell_counts.tsv',sep='\t');main=manifest[(manifest.config=='primary')&(manifest.cell_type=='Endothelial')].set_index('sample');samples=list(main.index)
 counts=pd.concat([pd.read_csv(folder/f'{s}_pseudobulk.tsv.gz',sep='\t',index_col=0)['primary__Endothelial'].rename(s) for s in samples],axis=1).fillna(0)
 logcpm=np.log2(counts/counts.sum(axis=0)*1e6+1)
 de=pd.read_csv(results/'primary__Endothelial__de.tsv.gz',sep='\t',index_col=0)
 for panel,spec in PANELS.items():
  available=[g for g in spec['genes'] if g in logcpm.index]
  sc=logcpm.loc[available].mean(axis=0)
  for sample in samples:scores.append(dict(dataset=acc,sample=sample,condition=main.loc[sample,'condition'],panel=panel,score=float(sc[sample]),unit='mean log2 pseudobulk CPM+1'))
  for gene in spec['genes']:
   if gene in de.index:
    r=de.loc[gene];bbb.append(dict(dataset=acc,panel=panel,gene=gene,effect=r.log2FoldChange,ci_low=r.log2FoldChange-1.96*r.lfcSE,ci_high=r.log2FoldChange+1.96*r.lfcSE,pvalue=r.pvalue,padj=r.padj,method='PyDESeq2 Wald',control='Sham'))
 for config in ['primary','strict_identity','mt10','low_myeloid']:
  table=pd.read_csv(results/f'{config}__communication.tsv.gz',sep='\t')
  sampletable=pd.read_csv(results/f'{config}__communication_sample_scores.tsv.gz',sep='\t')
  for sender,receiver,ligand,receptor in PAIRS:
   mask=(table.sender==sender)&(table.receiver==receiver)&(table.ligand==ligand)&(table.receptor==receptor)
   if not mask.any():continue
   r=table[mask].iloc[0].to_dict();st=sampletable[(sampletable.sender==sender)&(sampletable.receiver==receiver)&(sampletable.ligand==ligand)&(sampletable.receptor==receptor)].copy()
   loo=[]
   for s in st['sample']:
    sub=st[st['sample']!=s];loo.append(sub.loc[sub.condition=='MCAO','score'].mean()-sub.loc[sub.condition=='Sham','score'].mean())
   r.update(loso_min_delta=min(loo),loso_max_delta=max(loo),loso_same_direction_fraction=float((np.sign(loo)==np.sign(r['score_difference'])).mean()))
   pairrows.append(r)
   if config!='primary':continue
   for role,typ,complex_ in [('ligand',sender,ligand),('receptor',receiver,receptor)]:
    d=pd.read_csv(results/f'primary__{typ}__de.tsv.gz',sep='\t',index_col=0)
    for gene in complex_.split('_'):
     if gene in d.index:
      z=d.loc[gene];components.append(dict(dataset=acc,sender=sender,receiver=receiver,ligand=ligand,receptor=receptor,role=role,cell_type=typ,gene=gene,effect=z.log2FoldChange,ci_low=z.log2FoldChange-1.96*z.lfcSE,ci_high=z.log2FoldChange+1.96*z.lfcSE,pvalue=z.pvalue,padj=z.padj))
   for panel in ['BBB_junction','BBB_transport']:
    sc=pd.DataFrame(scores);sc=sc[(sc.dataset==acc)&(sc.panel==panel)]
    m=st.merge(sc,on=['sample','condition'],suffixes=('_communication','_panel'));xc=m.score_communication-m.groupby('condition').score_communication.transform('mean');yc=m.score_panel-m.groupby('condition').score_panel.transform('mean')
    associations.append(dict(dataset=acc,sender=sender,receiver=receiver,ligand=ligand,receptor=receptor,panel=panel,n_samples=len(m),spearman_unadjusted=float(stats.spearmanr(m.score_communication,m.score_panel).statistic),pearson_within_condition=float(stats.pearsonr(xc,yc).statistic),interpretation='descriptive only; low degrees of freedom; no causal inference'))
 # Sensitivity on the full common resource universe, no binomial p-values.
 baseline=pd.read_csv(results/'primary__communication.tsv.gz',sep='\t')
 for config in ['strict_identity','mt10','low_myeloid']:
  alt=pd.read_csv(results/f'{config}__communication.tsv.gz',sep='\t');m=baseline.merge(alt,on=['sender','receiver','ligand','receptor'],suffixes=('_primary','_alternative'));m=m[m.eligible_10_primary&m.eligible_10_alternative]
  sens.append(dict(dataset=acc,config=config,n_common=len(m),same_direction=int((np.sign(m.score_difference_primary)==np.sign(m.score_difference_alternative)).sum()),rho=float(stats.spearmanr(m.score_difference_primary,m.score_difference_alternative).statistic)))
sorted_de=pd.read_csv(ROOT/'sorted_rna/diff_expression.tsv.gz',sep='\t')
for panel,spec in PANELS.items():
 for gene in spec['genes']:
  z=sorted_de[(sorted_de.gene==gene)&(sorted_de.cell_type=='Endothelial')]
  if not len(z):continue
  r=z.iloc[0];bbb.append(dict(dataset='GSE163752',panel=panel,gene=gene,effect=r.log2_difference,ci_low=r.ci_low,ci_high=r.ci_high,pvalue=r.pvalue,padj=r.padj,method='paired t on log2 normalized+1',control='Contralateral'))
for sender,receiver,ligand,receptor in PAIRS:
 for role,typ,complex_ in [('ligand',sender,ligand),('receptor',receiver,receptor)]:
  for gene in complex_.split('_'):
   z=sorted_de[(sorted_de.gene==gene)&(sorted_de.cell_type==typ)]
   if not len(z):continue
   r=z.iloc[0];components.append(dict(dataset='GSE163752',sender=sender,receiver=receiver,ligand=ligand,receptor=receptor,role=role,cell_type=typ,gene=gene,effect=r.log2_difference,ci_low=r.ci_low,ci_high=r.ci_high,pvalue=r.pvalue,padj=r.padj))
for name,rows in [('barrier_gene_effects',bbb),('endothelial_program_samples',scores),('illustrative_pair_components',components),('illustrative_pair_stability',pairrows),('communication_barrier_associations',associations),('sensitivity_summary',sens)]:pd.DataFrame(rows).to_csv(OUT/f'{name}.tsv',sep='\t',index=False)
print(pd.DataFrame(bbb).pivot(index='gene',columns='dataset',values='effect').round(3).to_string());print(pd.DataFrame(sens).to_string(index=False))
