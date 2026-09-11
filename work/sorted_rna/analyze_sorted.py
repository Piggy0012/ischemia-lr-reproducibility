from pathlib import Path
import re,json
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
b=Path(__file__).resolve().parent
df=pd.read_csv(b/'GSE163752_counts.matrix.norm_anno.txt.gz',sep='\t',low_memory=False)
cols=[c for c in df if re.fullmatch('[cs](AC|EC)_[1-6]',c)]
assert len(cols)==24 and df[cols].notna().all().all() and (df[cols]>=0).all().all()
expr=df.dropna(subset=['Ensembl gene']).groupby('Ensembl gene')[cols].sum();expr.index.name='gene'
expr.to_csv(b/'astro_endothelial_normalized.tsv.gz',sep='\t')
manifest=[];res=[]
for typ,abbr in [('Astrocyte','AC'),('Endothelial','EC')]:
    control=[f'c{abbr}_{i}' for i in range(1,7)];stroke=[f's{abbr}_{i}' for i in range(1,7)]
    a=expr[control].to_numpy();c=expr[stroke].to_numpy();mask=(np.concatenate([a,c],axis=1)>=10).sum(axis=1)>=6
    genes=expr.index[mask];a=a[mask];c=c[mask];la=np.log2(a+1);lc=np.log2(c+1)
    paired=stats.ttest_rel(lc,la,axis=1);unpaired=stats.ttest_ind(lc,la,axis=1,equal_var=False)
    effect=(lc-la).mean(axis=1);se=(lc-la).std(axis=1,ddof=1)/np.sqrt(6)
    out=pd.DataFrame({'gene':genes,'cell_type':typ,'log2_difference':effect,'paired_t':paired.statistic,'pvalue':paired.pvalue,'welch_p':unpaired.pvalue,'ci_low':effect-stats.t.ppf(.975,5)*se,'ci_high':effect+stats.t.ppf(.975,5)*se,'control_mean':a.mean(axis=1),'stroke_mean':c.mean(axis=1),'n_pairs_same_direction':np.maximum(((lc-la)>0).sum(axis=1),((lc-la)<0).sum(axis=1))})
    out['padj']=multipletests(out.pvalue.fillna(1),method='fdr_bh')[1];out['welch_padj']=multipletests(out.welch_p.fillna(1),method='fdr_bh')[1]
    res.append(out)
    for condition,samples in [('Contralateral',control),('Ipsilateral',stroke)]:
        for s in samples:manifest.append({'sample':s,'cell_type':typ,'condition':condition,'replicate_index':int(s[-1]),'unit':'independent pool of 3-4 mice','pairing':'matched by deposited replicate index; unpaired sensitivity reported'})
pd.concat(res).to_csv(b/'diff_expression.tsv.gz',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(b/'sample_manifest.tsv',sep='\t',index=False)
audit={'rows_original':len(df),'genes_after_symbol_aggregation':len(expr),'samples_used':24,'cell_types_used':['Astrocyte','Endothelial'],'data_scale':'author-provided normalized counts; log2(x+1) for inference','model':'paired t-test on log2 normalized expression, BH within cell type; Welch sensitivity','control':'contralateral hemisphere, not sham','excluded_mapping_ambiguity':'12 ipsilateral pericyte/microglia GSM titles conflict with matrix descriptions; PC/MG excluded from formal validation','no_claim':'This analysis does not validate single-cell subtype identity, physical communication or causal barrier function.'}
(b/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(audit)
for out in res:print(out.cell_type.iloc[0],len(out),'FDR05',int((out.padj<.05).sum()))
