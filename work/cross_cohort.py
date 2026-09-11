from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/cross_cohort';OUT.mkdir(exist_ok=True)
stats={}
for typ in ['Astrocyte','Endothelial']:
 a=pd.read_csv(ROOT/f'results/GSE174574/primary__{typ}__de.tsv.gz',sep='\t',index_col=0)
 b=pd.read_csv(ROOT/f'results/GSE245386/primary__{typ}__de.tsv.gz',sep='\t',index_col=0)
 c=a.join(b,lsuffix='_discovery',rsuffix='_validation',how='inner')
 c['same_direction']=np.sign(c.log2FoldChange_discovery)==np.sign(c.log2FoldChange_validation)
 c['discovery_selected']=(c.padj_discovery<.05)&(c.log2FoldChange_discovery.abs()>=.5)
 c['both_selected']=c.discovery_selected&(c.padj_validation<.05)&(c.log2FoldChange_validation.abs()>=.5)
 c.to_csv(OUT/f'{typ}__gene_replication.tsv.gz',sep='\t')
 s=c[c.discovery_selected];both=c[c.both_selected & c.same_direction]
 stats[typ]={'common_tested':len(c),'rho_all':float(spearmanr(c.log2FoldChange_discovery,c.log2FoldChange_validation,nan_policy='omit').statistic),'discovery_selected_common':len(s),'same_direction_selected':int(s.same_direction.sum()),'same_direction_both_fdr_effect':len(both)}
 print(typ,stats[typ]);print(both.sort_values('padj_discovery')[['log2FoldChange_discovery','log2FoldChange_validation','padj_validation']].head(15).to_string())
keys=['sender','receiver','ligand','receptor']
a=pd.read_csv(ROOT/'results/GSE174574/primary__communication.tsv.gz',sep='\t');b=pd.read_csv(ROOT/'results/GSE245386/primary__communication.tsv.gz',sep='\t')
c=a.merge(b,on=keys,suffixes=('_discovery','_validation'),how='outer');c['same_direction']=np.sign(c.score_difference_discovery)==np.sign(c.score_difference_validation)
c['both_eligible10']=c.eligible_10_discovery.fillna(False)&c.eligible_10_validation.fillna(False)
c['discovery_selected']=c.eligible_10_discovery.fillna(False)&c.candidate_changed_component_discovery.fillna(False)
c.to_csv(OUT/'communication_replication.tsv.gz',sep='\t',index=False)
s=c[c.both_eligible10];d=s[s.discovery_selected]
stats['communication']={'discovery_eligible10':int(a.eligible_10.sum()),'validation_eligible10':int(b.eligible_10.sum()),'common_eligible10':len(s),'common_same_direction':int(s.same_direction.sum()),'discovery_selected':int((a.eligible_10&a.candidate_changed_component).sum()),'selected_common':len(d),'selected_common_same_direction':int(d.same_direction.sum()),'rho_all_common':float(spearmanr(s.score_difference_discovery,s.score_difference_validation).statistic)}
print(stats['communication'])
print(d[d.same_direction].sort_values('score_difference_discovery',key=lambda x:x.abs(),ascending=False)[keys+['score_difference_discovery','score_difference_validation','welch_fdr_eligible_10_validation','ligand_mean_log2fc_validation','receptor_mean_log2fc_validation']].head(30).to_string(index=False))
(OUT/'summary.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
