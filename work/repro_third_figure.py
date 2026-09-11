"""Figure 5. Preview and actual visual inspection must precede --export."""
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent/'outputs/reproducibility_v2'
FIG = OUT/'figures'
DATA = OUT/'figure_source_data'
BASE = 'Figure5_third_cohort_transfer'
METRICS = ['custom_score','native_magnitude_priority','diagnostic_unique_column_priority']
LABELS = ['Custom','Native RRA','Unique-column RRA']
SKILL = Path.home()/'.codex/skills/scipilot-figure-skill/scripts'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':400})

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def build(export=False):
    for p in [FIG,DATA]: p.mkdir(parents=True,exist_ok=True)
    qp = ROOT/'repro_third_cohort/library_qc_summary.tsv'
    rp = OUT/'tables/third_rank_diagnostic_summary.tsv'
    q = pd.read_csv(qp,sep='\t',float_precision='round_trip')
    full = pd.read_csv(rp,sep='\t',float_precision='round_trip')
    d = full[(full.universe=='all_three_cohorts_complete') & full.metric.isin(METRICS)].copy()
    assert len(q)==6 and q.condition.value_counts().to_dict()=={'MCAO':3,'Sham':3}
    assert len(d)==6 and not d.duplicated(['baseline_dataset','metric']).any()
    assert d.n_candidates.nunique()==1 and d.n_candidates.iloc[0]>0
    assert set(d.config)=={'primary'}
    assert np.allclose(d.n_same_direction/d.n_candidates,d.same_direction_all_fraction)
    d['n_zero_in_either'] = d.n_candidates-d.n_nonzero_both
    count_data=q[['sample','condition','library_replicate','n_qc','Astrocyte','Endothelial',
                  'n_qc_unassigned','qc_unassigned_fraction']].copy()
    count_data.to_csv(DATA/'figure5_library_counts.tsv',sep='\t',index=False)
    d.to_csv(DATA/'figure5_fixed_17_library_comparison.tsv',sep='\t',index=False)
    # Preserve every method and pairwise/global universe in the same source-data folder.
    full.to_csv(DATA/'figure5_all_metric_and_universe_results.tsv',sep='\t',index=False)
    profile={'unit_counts':'one original-primary assigned-nucleus count per deposited library',
        'unit_effects':'aggregate cross-cohort comparisons; candidate rows share genes and are not animal replicates',
        'biological_replication':'Table S1 reports 3 biological replicates and 3 libraries/group; mice per library unspecified',
        'choice':'individual library dot plots; unconnected rank-correlation dots with explicit direction/zero denominators',
        'alternative':'annotated numerical table; less immediate visual comparison of correlation',
        'no_new_pvalues_or_uncertainty_intervals':True,
        'shared_universe':'same native/diagnostic complete candidates across all 17 libraries and all metrics',
        'third_count_source':'raw only','third_selection':'original primary only',
        'input_sha256':{str(qp.relative_to(ROOT)):sha(qp),str(rp.relative_to(ROOT.parent)):sha(rp)}}
    if SKILL.exists():
        sys.path.insert(0,str(SKILL))
        from profile_data import profile_data
        profile['counts_profile']=profile_data(str(DATA/'figure5_library_counts.tsv'),group_cols=['condition'])
        profile['comparison_profile']=profile_data(str(DATA/'figure5_fixed_17_library_comparison.tsv'),group_cols=['baseline_dataset'])
    (DATA/'figure5_data_profile.json').write_text(json.dumps(profile,indent=2,default=str),encoding='utf-8')
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.3),layout='constrained',gridspec_kw={'height_ratios':[1,1.15]})
    colors={'Sham':'#0072B2','MCAO':'#D55E00'}
    markers={'Sham':'o','MCAO':'^'}
    for ax,typ,title,limit in zip(axes[0],['Astrocyte','Endothelial'],['Assigned astrocyte nuclei','Assigned endothelial nuclei'],[1500,300]):
        for j,condition in enumerate(['Sham','MCAO']):
            part=q[q.condition==condition].sort_values('library_replicate')
            xs=j+np.array([-.14,0,.14])
            for x,row in zip(xs,part.itertuples()):
                y=getattr(row,typ)
                ax.scatter(x,y,color=colors[condition],marker=markers[condition],s=25,zorder=3)
                ax.annotate(('S' if condition=='Sham' else 'M')+str(row.library_replicate),(x,y),
                            xytext=(0,5),textcoords='offset points',ha='center',va='bottom',fontsize=7)
        ax.set_xlim(-.5,1.5); ax.set_ylim(0,limit)
        ax.set_xticks([0,1],['Sham\n3 libraries','MCAO 24 h\n3 libraries'])
        ax.set_ylabel('Nuclei per library'); ax.set_title(title)
        ax.grid(axis='y',color='#E7E7E7',lw=.45,zorder=0)
        ax.set_axisbelow(True)
    n=int(d.n_candidates.iloc[0])
    finite=d.spearman_rho.dropna()
    minimum=-1 if (finite<0).any() else 0
    for ax,baseline in zip(axes[1],['GSE174574','GSE245386']):
        b=d[d.baseline_dataset==baseline].set_index('metric').loc[METRICS]
        labels=[]
        for j,(label,row) in enumerate(zip(LABELS,b.itertuples())):
            labels.append(f'{label}\n{row.n_same_direction}/{n} agree; {row.n_zero_in_either} zero')
            if pd.notna(row.spearman_rho):
                ax.scatter(row.spearman_rho,j,s=29,marker=['o','s','D'][j],color='#242424',zorder=3)
                off=.045 if row.spearman_rho<.85 else -.045
                ax.text(row.spearman_rho+off,j,f'{row.spearman_rho:.3f}',va='center',
                        ha='left' if off>0 else 'right',fontsize=7)
            else:
                ax.text(.5,j,'undefined',va='center',ha='center',fontsize=7)
        ax.set_xlim(minimum,1); ax.set_ylim(2.55,-.55)
        ax.set_yticks(range(3),labels)
        ax.set_xticks([0,.25,.5,.75,1] if minimum==0 else [-1,-.5,0,.5,1])
        ax.set_xlabel('Cross-cohort Spearman rho')
        ax.set_title('GSE332910 vs '+baseline)
        ax.grid(axis='x',color='#E7E7E7',lw=.45,zorder=0)
    fig.suptitle(f'Acute striatal nuclear extension | fixed N = {n} candidates across 17 libraries',fontsize=9)
    fig.canvas.draw()
    if SKILL.exists():
        from layout_tools import add_panel_labels
        add_panel_labels(fig,style='upper',fontsize=9)
    else:
        for label,ax in zip('ABCD',axes.flat):
            p=ax.get_position();fig.text(p.x0-.04,p.y1+.025,label,fontweight='bold',fontsize=9)
    fig.canvas.draw()
    issues=[]
    if SKILL.exists():
        from visual_qa import audit_layout
        issues=audit_layout(fig)
    (DATA/'figure5_layout_audit.json').write_text(json.dumps(issues,indent=2,default=str),encoding='utf-8')
    png=FIG/(BASE+'.png')
    fig.savefig(png,dpi=400)
    Image.open(png).convert('L').save(DATA/'figure5_grayscale_preview.png',dpi=(400,400))
    if export:
        review_path=DATA/'figure5_visual_review.json'
        review=json.loads(review_path.read_text(encoding='utf-8'))
        assert review['passed'] is True and review['preview_sha256']==sha(png), 'Inspect the current preview before vector export'
        fig.savefig(FIG/(BASE+'.pdf'))
        fig.savefig(FIG/(BASE+'.svg'))
    plt.close(fig)
    print(json.dumps({'preview':str(png),'sha256':sha(png),'N':n,'exported_vectors':export,'issues':issues},indent=2,default=str))

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--export',action='store_true')
    build(ap.parse_args().export)
