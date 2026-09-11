"""Descriptive paired-pool effects; no new significance classification."""
from pathlib import Path
import json,sys,hashlib
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2';TABLES=OUT/'tables';FIG=OUT/'figures';DATA=OUT/'figure_source_data'
for p in [TABLES,FIG,DATA]:p.mkdir(parents=True,exist_ok=True)
source=ROOT/'sorted_rna/astro_endothelial_normalized.tsv.gz'
expr=pd.read_csv(source,sep='\t',index_col=0)
old=pd.read_csv(ROOT/'sorted_rna/diff_expression.tsv.gz',sep='\t')
selected={'Astrocyte':['Spp1','Timp3','Ptprz1','Lrp1','Col4a1'],
          'Endothelial':['Itga5','Itgb1','Kdr','Ptn','Plat','Itga3']}
summary=[];points=[];maxerror=0
for typ,abbr in [('Astrocyte','AC'),('Endothelial','EC')]:
    a=expr[[f'c{abbr}_{i}' for i in range(1,7)]].to_numpy()
    b=expr[[f's{abbr}_{i}' for i in range(1,7)]].to_numpy()
    diff=np.log2(b+1)-np.log2(a+1)
    mean=diff.mean(axis=1);se=diff.std(axis=1,ddof=1)/np.sqrt(6)
    ci=stats.t.ppf(.975,5)*se
    tested=(np.concatenate([a,b],axis=1)>=10).sum(axis=1)>=6
    t=pd.DataFrame({'gene':expr.index,'cell_type':typ,'paired_log2_effect':mean,
        'ci_low':mean-ci,'ci_high':mean+ci,'n_paired_pools':6,'previous_transcriptome_filter':tested,
        'n_positive_pairs':(diff>0).sum(axis=1),'n_negative_pairs':(diff<0).sum(axis=1),
        'interpretation':'Author-normalized expression; paired ipsilateral-minus-contralateral log2(x+1); pointwise t interval'})
    o=old[old.cell_type==typ].set_index('gene')
    j=t.set_index('gene').loc[o.index]
    error=max(np.abs(j.paired_log2_effect-o.log2_difference).max(),np.abs(j.ci_low-o.ci_low).max(),np.abs(j.ci_high-o.ci_high).max())
    maxerror=max(maxerror,float(error));assert error<1e-10
    summary.append(t)
    for g in selected[typ]:
        assert g in expr.index
        i=expr.index.get_loc(g)
        for k in range(6):points.append({'cell_type':typ,'gene':g,'pool_pair':k+1,'paired_log2_difference':diff[i,k],
            'control_author_normalized':a[i,k],'stroke_author_normalized':b[i,k]})
full=pd.concat(summary,ignore_index=True)
full.to_csv(TABLES/'sorted_all_gene_descriptive_effects.tsv.gz',sep='\t',index=False)
sel=pd.concat([full[(full.cell_type==typ)&full.gene.isin(gs)] for typ,gs in selected.items()],ignore_index=True)
sel.to_csv(TABLES/'sorted_fixed_component_effects_ci.tsv',sep='\t',index=False)
pd.DataFrame(points).to_csv(DATA/'sorted_fixed_components_per_pool.tsv',sep='\t',index=False)
audit={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'old_effect_and_ci_max_absolute_error':maxerror,'n_paired_pools':6,'pool_size_mice':'3-4',
    'paired_by':'Deposited replicate number','control':'contralateral hemisphere',
    'new_p_values_computed':False,'significance_classification':False,
    'all_original_tests_retained_at':'work/sorted_rna/diff_expression.tsv.gz',
    'ci':'Pointwise two-sided 95% t interval, df=5; not simultaneous, normality of paired differences assumed',
    'selection':'Five illustrative LR examples fixed in prior revision; descriptive illustration, not prespecified hypotheses'}
(TABLES/'sorted_descriptive_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'xtick.labelsize':7,
    'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(7.2,3.8),layout='constrained')
pointdf=pd.DataFrame(points)
for k,(typ,gs) in enumerate(selected.items()):
    ax=axes[k]
    for i,g in enumerate(gs):
        r=sel[(sel.cell_type==typ)&(sel.gene==g)].iloc[0]
        d=pointdf[(pointdf.cell_type==typ)&(pointdf.gene==g)]
        ax.scatter(d.paired_log2_difference,i+np.linspace(-.14,.14,6),s=16,color='#0072B2',alpha=.7,linewidths=0)
        ax.errorbar(r.paired_log2_effect,i,xerr=[[r.paired_log2_effect-r.ci_low],[r.ci_high-r.paired_log2_effect]],
                    color='black',fmt='s',ms=3,capsize=2,lw=1)
    ax.axvline(0,color='#888888',ls='--',lw=.7);ax.set_yticks(range(len(gs)),gs);ax.invert_yaxis()
    ax.set_title(typ+'\nSix paired pools; 3-4 mice per pool')
    ax.set_xlabel('Ipsilateral - contralateral log2(x+1)')
    ax.text(-.12,1.04,'AB'[k],transform=ax.transAxes,fontweight='bold',fontsize=11)
fig.savefig(FIG/'FigureS4_sorted_effects_and_intervals.png',dpi=400)
if '--export' in sys.argv:
    for ext in ['pdf','svg']:fig.savefig(FIG/f'FigureS4_sorted_effects_and_intervals.{ext}')
print(json.dumps({'n_descriptive_rows':len(full),'components':len(sel),'prior_effect_ci_error':maxerror}))
