"""Figure 3: preview/inspect PNG first; pass --export only after visual review."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent/'outputs/reproducibility_v2'
FIG = OUT/'figures'; DATA = OUT/'figure_source_data'
for p in [FIG, DATA]: p.mkdir(parents=True, exist_ok=True)
BASE = 'Figure3_implementation_and_rank_universe'
EXPORT = '--export' in sys.argv
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,
    'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,
    'axes.spines.right':False,'savefig.dpi':400})
BLUE='#0072B2'; ORANGE='#D55E00'; GRAY='#B6B6B6'
METRICS=['native_full_network_priority','unique_full_network_priority','unique_fixed_global_priority']
LABELS=['Native v1.10.0','Unique-column RRA','+ fixed global set']
SKILL=Path.home()/'.codex/skills/scipilot-figure-skill/scripts'

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def build():
    path=OUT/'tables/rank_background_summary.tsv'
    d=pd.read_csv(path,sep='\t')
    assert len(d)==6 and not d.duplicated(['config','metric']).any() and d.notna().all().all()
    assert set(d.metric)==set(METRICS)
    assert np.allclose(d.n_same_direction/d.n_candidates,d.same_direction_all_fraction)
    assert np.allclose(d.n_same_direction/d.n_nonzero_both,d.same_direction_nonzero_fraction)
    d['n_zero_in_either_cohort']=d.n_candidates-d.n_nonzero_both
    assert (d.n_zero_in_either_cohort>0).all()
    context=pd.read_csv(OUT/'tables/rank_background_context.tsv',sep='\t')
    d.to_csv(DATA/'figure3_comparison_and_denominators.tsv',sep='\t',index=False)
    context.to_csv(DATA/'figure3_global_universe_context.tsv',sep='\t',index=False)
    pd.read_csv(OUT/'tables/rank_background_zero_effects.tsv',sep='\t').to_csv(DATA/'figure3_zero_effect_sources.tsv',sep='\t',index=False)
    profile={'shape':list(d.shape),'missing_values':int(d.isna().sum().sum()),
             'unit':'one aggregate comparison per implementation/selection; rows are not animal replicates',
             'chart_choice':'unconnected horizontal dots and composition bars',
             'alternative':'annotated table or heatmap; less direct for showing denominator loss',
             'selection_count':2,'fixed_target_counts':[187,174],
             'point_estimates_have_no_error_bars_or_new_pvalues':True}
    if SKILL.exists():
        sys.path.insert(0,str(SKILL))
        from profile_data import profile_data
        profile['skill_profile']=profile_data(str(path),group_cols=['config'])
    (DATA/'figure3_data_profile.json').write_text(json.dumps(profile,indent=2,default=str),encoding='utf-8')
    fig,axes=plt.subplots(3,2,figsize=(7.2,7.2),layout='constrained',
                          gridspec_kw={'height_ratios':[.92,1.15,1.0]})
    for col,config in enumerate(['primary','reference_singlet']):
        b=d[d.config==config].set_index('metric').loc[METRICS]
        color,marker=(BLUE,'o') if col==0 else (ORANGE,'^')
        n=int(b.n_candidates.iloc[0])
        ax=axes[0,col]
        for j,r in enumerate(b.itertuples()):
            ax.scatter(r.spearman_rho,j,color=color,marker=marker,s=30,zorder=3)
            ax.text(r.spearman_rho+.04,j,f'{r.spearman_rho:.3f}',va='center',fontsize=7)
        ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.set_ylim(2.5,-.5)
        ax.set_yticks(range(3),LABELS);ax.set_xlabel('Cross-cohort Spearman rho')
        ax.set_title(('Primary selection' if col==0 else 'Reference-supported singlets')+f' | N = {n}')
        ax.grid(axis='x',color='#E7E7E7',lw=.45,zorder=0)
        ax=axes[1,col]
        for j,r in enumerate(b.itertuples()):
            for off,x,fill in [(-.13,r.same_direction_all_fraction*100,True),
                                (.13,r.same_direction_nonzero_fraction*100,False)]:
                ax.scatter(x,j+off,facecolor=color if fill else 'white',edgecolor=color,
                           linewidth=1,marker=marker,s=29,zorder=3)
                tx=x+3 if x<76 else x-3
                ax.text(tx,j+off,f'{x:.1f}',ha='left' if x<76 else 'right',va='center',fontsize=7)
        ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.set_ylim(2.95,-.5)
        ax.set_yticks(range(3),LABELS);ax.set_xlabel('Same-direction effects (%)')
        ax.set_title('All candidates and nonzero candidates')
        ax.grid(axis='x',color='#E7E7E7',lw=.45,zorder=0)
        ax.legend(handles=[Line2D([],[],lw=0,marker=marker,markersize=4,markeredgecolor=color,
                                  markerfacecolor=color,label='All N'),
                           Line2D([],[],lw=0,marker=marker,markersize=4,markeredgecolor=color,
                                  markerfacecolor='white',label='Nonzero in both')],
                  loc='lower left',frameon=False,ncol=2,handletextpad=.35,columnspacing=.8,
                  borderpad=.15,borderaxespad=.05)
        ax=axes[2,col]
        for j,r in enumerate(b.itertuples()):
            good=100*r.n_nonzero_both/r.n_candidates
            zero=100-good
            ax.barh(j,good,color=color,alpha=.28,height=.5,edgecolor='white',linewidth=.5)
            ax.barh(j,zero,left=good,color=GRAY,height=.5,edgecolor='white',linewidth=.5,hatch='///')
            ax.text(good/2,j,str(int(r.n_nonzero_both)),ha='center',va='center',fontsize=7)
            ax.text(good+zero/2,j,str(int(r.n_zero_in_either_cohort)),ha='center',va='center',fontsize=7)
        ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.set_ylim(2.95,-.5)
        ax.set_yticks(range(3),LABELS);ax.set_xlabel('Share of the same N candidates (%)')
        ax.set_title('Nonzero disease effects; counts shown')
        ax.legend(handles=[Patch(facecolor=color,alpha=.28,label='Nonzero in both'),
                           Patch(facecolor=GRAY,hatch='///',label='Zero in either')],
                  loc='lower left',frameon=False,ncol=2,handlelength=1,handletextpad=.4,
                  columnspacing=.8,borderpad=.15,borderaxespad=.05)
    for i,ax in enumerate(axes.flat):
        ax.text(-.15,1.06,'ABCDEF'[i],transform=ax.transAxes,fontweight='bold',fontsize=10)
    fig.get_layout_engine().set(h_pad=.12,w_pad=.04,hspace=.07,wspace=.08)
    issues=[]
    if SKILL.exists():
        from visual_qa import audit_layout
        issues=audit_layout(fig)
    if any(sev=='FAIL' for sev,_ in issues):
        raise AssertionError(issues)
    png=FIG/(BASE+'.png')
    fig.savefig(png,dpi=400)
    gray=FIG/(BASE+'_grayscale.png')
    with Image.open(png) as im: im.convert('L').save(gray,dpi=(400,400))
    outputs={'png':digest(png),'grayscale_png':digest(gray)}
    if EXPORT:
        for ext in ['pdf','svg']:
            p=FIG/(BASE+'.'+ext);fig.savefig(p);outputs[ext]=digest(p)
    qa={'status':'exported_after_review' if EXPORT else 'png_ready_for_visual_review',
        'size_inches':[7.2,7.2],'dpi':400,'font':'Arial','minimum_font_size_pt':7,
        'semantic_review':{'no_new_pvalues':True,'no_error_bars':True,'both_direction_denominators_shown':True,
                           'same_fixed_target_n_across_stages':True,'no_lines_connecting_categorical_scores':True,
                           'composition_bars_are_edge_counts_not_animal_means':True},
        'machine_layout_issues':issues,'input_sha256':{path.name:digest(path)},
        'script_sha256':digest(Path(__file__)),'output_sha256':outputs}
    (DATA/'figure3_qa.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False),encoding='utf-8')
    plt.close(fig)
    print(json.dumps(qa,ensure_ascii=False))

if __name__=='__main__':build()
