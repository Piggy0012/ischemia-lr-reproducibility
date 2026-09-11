"""S6 controlled CellChat extension: inspect preview before --export vectors."""
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from PIL import Image

ROOT=Path(__file__).resolve().parent
TABLE=ROOT/'repro_v3_cellchat/tables'
OUT=ROOT.parent/'outputs/reproducibility_v3_cellchat'
FIG=OUT/'figures';DATA=OUT/'figure_source_data'
BASE='FigureS6_cellchat_controlled_extension'
KEY=['source','target','ligand','receptor']
METRICS=['custom_coavailability','cellchat_probability','cellchat_strength_percentile','LIANA_unique_column_RRA_priority']
TITLES=['Expression coavailability','CellChat native inferred strength (net$prob)',
        'CellChat strength percentile (analyst-derived)','LIANA unique-column RRA priority']
SKILL=Path.home()/'.codex/skills/scipilot-figure-skill/scripts'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
  'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,
  'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':400})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
def write_json(p,d):p.write_text(json.dumps(d,indent=2,default=str),encoding='utf-8')

def build(export=False):
    for p in [FIG,DATA]:p.mkdir(parents=True,exist_ok=True)
    paths=[TABLE/n for n in ['target_disease_effects.tsv.gz','cross_cohort_concordance.tsv',
      'zero_strength_target_audit.tsv','tie_zero_audit.tsv','summary_audit.json']]
    d=read(paths[0]);summ=read(paths[1]);zero=read(paths[2])
    z=zero.set_index(KEY).all_11_strength_zero
    s=summ[summ.scope=='all_fixed_targets'].set_index('metric').loc[METRICS]
    assert (s.n_candidates==32).all() and int(z.sum())==3
    points=[]
    for metric in METRICS:
        p=d[d.metric==metric].pivot(index=KEY,columns='dataset',values='mcao_minus_sham')
        assert len(p)==32 and np.isfinite(p.to_numpy()).all()
        p['all_11_native_strength_zero']=z.loc[p.index]
        p['metric']=metric;points.append(p.reset_index())
    source=pd.concat(points,ignore_index=True)
    source.to_csv(DATA/'figureS6_candidate_effects.tsv',sep='\t',index=False)
    s.reset_index().to_csv(DATA/'figureS6_concordance.tsv',sep='\t',index=False)
    read(paths[3]).to_csv(DATA/'figureS6_zero_tie_floors.tsv',sep='\t',index=False)
    profile={'question':'Does expression-effect concordance transfer to native CellChat strength and analyst-derived relative priority in the exactly matched subset?',
      'observation_unit':'one directed ligand-receptor candidate, not an independent animal replicate',
      'fixed_target_N':32,'fixed_network_N':1222,'input_libraries':11,
      'resource_scope':'32 of 187 original targets have exact canonical native mouse resource definitions; the remaining 155 are outside this controlled subset',
      'recommended_chart':'Four scatter panels retain all candidate effects, zeros, direction classes and shared all-zero-strength flags.',
      'alternatives':['A compact metrics table is equally exact but conceals ties and zero effects.',
                      'A rho dot plot is compact but conceals the native zero-effect coverage difference.'],
      'scales':'Linear raw group differences; each axis padded to its observed range with zero included. Different metrics have different units.',
      'overlap':'No jitter is applied. Identical candidate coordinates may overlap, including three all-zero-strength candidates.',
      'no_fitted_line_or_inference':'Candidate dependencies and small animal groups preclude treating plotted points as independent replicates; no P or CI displayed.',
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    # Optional local authoring helpers are not a public runtime dependency.
    profile_data=None; audit_layout=None
    if SKILL.exists():
        sys.path.insert(0,str(SKILL))
        try:
            from profile_data import profile_data
            from visual_qa import audit_layout
        except ImportError:
            profile_data=None; audit_layout=None
    profile['optional_authoring_helper_unavailable']=profile_data is None or audit_layout is None
    profile['data_profile']=profile_data(str(DATA/'figureS6_candidate_effects.tsv'),group_cols=['metric','source']) if profile_data else {
      'rows':len(source),'missing_values':source.isna().sum().to_dict(),
      'numeric_summary':source[['GSE174574','GSE245386']].describe().to_dict()}
    write_json(DATA/'figureS6_data_profile.json',profile)

    fig,axes=plt.subplots(2,2,figsize=(7.2,6.2))
    fig.subplots_adjust(left=.095,right=.98,bottom=.115,top=.83,wspace=.33,hspace=.61)
    colors={'Astrocyte':'#0072B2','Endothelial':'#D55E00'}
    markers={'Astrocyte':'o','Endothelial':'^'}
    for letter,ax,metric,title in zip('ABCD',axes.flat,METRICS,TITLES):
        p=source[source.metric==metric]
        for hollow in [False,True]:
            for direction in colors:
                g=p[(p.source==direction)&(p.all_11_native_strength_zero==hollow)]
                ax.scatter(g.GSE174574,g.GSE245386,s=23 if not hollow else 36,
                  marker=markers[direction],edgecolor=colors[direction],linewidth=.8,
                  facecolor='none' if hollow else colors[direction],alpha=1 if hollow else .82,
                  zorder=4 if hollow else 3)
        for dim,col in [('x','GSE174574'),('y','GSE245386')]:
            low=min(float(p[col].min()),0);high=max(float(p[col].max()),0)
            span=high-low if high>low else 1
            getattr(ax,'set_'+dim+'lim')(low-.10*span,high+.14*span)
            getattr(ax,dim+'axis').set_major_locator(MaxNLocator(nbins=5))
        ax.axhline(0,color='#AAAAAA',lw=.55,ls='--',zorder=1)
        ax.axvline(0,color='#AAAAAA',lw=.55,ls='--',zorder=1)
        row=s.loc[metric];same=int(row.n_same_direction);nz=int(row.n_nonzero_both)
        ax.set_title(title,y=1.12,pad=2)
        ax.text(0,1.045,f'ρ = {row.spearman_rho:.3f}   Same: {same}/32; nonzero: {same}/{nz}   Zero: {32-nz}',
                transform=ax.transAxes,ha='left',fontsize=7)
        ax.text(-.15,1.14,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
        ax.set_xlabel('GSE174574: mean MCAO − mean Sham')
        ax.set_ylabel('GSE245386: mean MCAO − mean Sham')
        if metric=='cellchat_strength_percentile':
            zz=p[p.all_11_native_strength_zero]
            assert zz.GSE174574.nunique()==zz.GSE245386.nunique()==1
            xx=float(zz.GSE174574.iloc[0]);yy=float(zz.GSE245386.iloc[0])
            ax.annotate('3 all-zero targets\n(overlap)',xy=(xx,yy),xytext=(.64,.83),
                        textcoords='axes fraction',fontsize=7,ha='center',va='bottom',
                        arrowprops={'arrowstyle':'-','color':'#777777','lw':.6},zorder=5)
    handles=[Line2D([],[],marker=markers[k],color=colors[k],linestyle='none',markersize=4.5,
                    label=label) for k,label in [('Astrocyte','Astrocyte → endothelial (8)'),('Endothelial','Endothelial → astrocyte (24)')]]
    fig.suptitle('Controlled shared resource: 32 targets; 1,222-edge rank background',y=.992,fontsize=9)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.968),ncol=2,frameon=False)
    fig.text(.5,.037,'Hollow symbols: native CellChat strength = 0 in all 11 libraries (3 targets).',ha='center',fontsize=7)
    fig.text(.5,.012,'Same = concordant nonzero effects; nonzero denominator = nonzero in both cohorts; Zero = zero in either.',ha='center',fontsize=7)
    fig.canvas.draw();issues=audit_layout(fig) if audit_layout else []
    write_json(DATA/'figureS6_layout_audit.json',issues)
    png=FIG/(BASE+'.png');fig.savefig(png,dpi=400)
    Image.open(png).convert('L').save(DATA/'figureS6_grayscale_preview.png',dpi=(400,400))
    if export:
        review=json.loads((DATA/'figureS6_visual_review.json').read_text(encoding='utf-8'))
        assert review['passed'] and review['preview_sha256']==sha(png),'Inspect this PNG before vector export'
        fig.savefig(FIG/(BASE+'.pdf'));fig.savefig(FIG/(BASE+'.svg'))
    plt.close(fig)
    print(json.dumps({'preview':str(png),'sha256':sha(png),'issues':issues,'vectors_exported':export},indent=2,default=str))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--export',action='store_true');build(ap.parse_args().export)
