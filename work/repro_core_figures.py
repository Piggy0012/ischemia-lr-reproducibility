"""Core reproducibility figures; preview first, --export only after visual review."""
from pathlib import Path
import argparse, json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent/'outputs/reproducibility_v2'
FIG = OUT/'figures'; DATA = OUT/'figure_source_data'
for p in [FIG, DATA]: p.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,
    'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,
    'axes.spines.right':False,'savefig.dpi':400})
BLUE='#0072B2'; ORANGE='#D55E00'; GRAY='#888888'; PURPLE='#CC79A7'
EXPORT = '--export' in sys.argv

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=400)
    if EXPORT:
        for ext in ['pdf','svg']: fig.savefig(FIG/(name+'.'+ext))
    plt.close(fig)

def label(ax,s): ax.text(-.12,1.06,s,transform=ax.transAxes,fontweight='bold',fontsize=11)

def figure1():
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.65),layout='constrained',
                           gridspec_kw={'height_ratios':[1,1.12]})
    rows=[]
    markers=['Slc1a2','Slc1a3','Aldh1l1','Pecam1','Cdh5','Kdr','Ptprc','Tyrobp','C1qa','Lyz2']
    audit=pd.read_csv(ROOT/'results/identity_marker_audit.tsv',sep='\t')
    for k,acc in enumerate(['GSE174574','GSE245386']):
        c=pd.read_csv(ROOT/'processed'/acc/'sample_cell_counts.tsv',sep='\t')
        c=c[(c.config=='primary')&c.cell_type.isin(['Astrocyte','Endothelial'])].copy()
        c['dataset']=acc; rows.append(c)
        ax=axes[0,k]
        for j,typ in enumerate(['Astrocyte','Endothelial']):
            for ci,cond in enumerate(['Sham','MCAO']):
                y=c.loc[(c.cell_type==typ)&(c.condition==cond),'n_cells'].to_numpy()
                xpos=3*j+ci
                ax.scatter(xpos+np.linspace(-.11,.11,len(y)),y,
                           c=BLUE if ci==0 else ORANGE,marker='o' if ci==0 else '^',
                           s=26,edgecolor='white',linewidth=.4,zorder=3)
                ax.plot([xpos-.17,xpos+.17],[y.mean()]*2,color='black',lw=.8)
        ax.set_xticks([0,1,3,4],['Sham','MCAO','Sham','MCAO']);ax.set_ylim(bottom=0)
        ax.set_ylabel('Cells per animal');ax.set_title(acc+('\n3 Sham / 3 MCAO' if k==0 else '\n3 Sham / 2 MCAO'))
        ax.text(.15,-.23,'Astrocytes',transform=ax.transAxes)
        ax.text(.64,-.23,'Endothelium',transform=ax.transAxes);label(ax,'AB'[k])
        d=audit[(audit.dataset==acc)&(audit.config=='primary')]
        t=d.groupby(['cell_type','gene']).fraction.mean().unstack('gene').reindex(
            index=['Astrocyte','Endothelial'],columns=markers)
        assert t.shape==(2,10) and t.notna().all().all()
        ax=axes[1,k]
        im=ax.imshow(t,vmin=0,vmax=1,cmap='cividis',aspect='auto')
        ax.axvline(5.5,color='white',lw=2)
        ax.set_xticks(range(10),markers,rotation=60,ha='right')
        ax.set_yticks([0,1],['Astrocyte','Endothelial'])
        ax.text(2.5,1.03,'Target lineage',ha='center',transform=ax.get_xaxis_transform(),fontsize=8)
        ax.text(7.5,1.03,'Myeloid',ha='center',transform=ax.get_xaxis_transform(),fontsize=8)
        label(ax,'CD'[k]);fig.colorbar(im,ax=ax,shrink=.68,label='Mean detection fraction')
        t.to_csv(DATA/f'figure1_markers_{acc}.tsv',sep='\t')
    pd.concat(rows).to_csv(DATA/'figure1_animal_counts.tsv',sep='\t',index=False)
    save(fig,'Figure1_samples_and_marker_groups')

def figure2():
    tab=OUT/'tables'
    effects=pd.read_csv(tab/'null_fixed_common_observed_effects.tsv.gz',sep='\t')
    summary=pd.read_csv(tab/'null_fixed_common_summary.tsv',sep='\t')
    dyn=pd.read_csv(tab/'null_dynamic_gate_all_allocations.tsv.gz',sep='\t')
    ds=pd.read_csv(tab/'null_dynamic_gate_summary.tsv',sep='\t')
    fig,axes=plt.subplots(2,2,figsize=(7.2,6.8),layout='constrained')
    # The actual column names are validated before plotting; no derived P values.
    for k,metric in enumerate(['custom_score','magnitude_priority']):
        long=effects[(effects.config=='primary')&(effects.metric==metric)]
        d=long.pivot(index=['sender','receiver','ligand','receptor'],columns='dataset',values='observed_mcao_minus_sham')
        assert len(d)==187
        ax=axes[0,k]
        x=d.GSE174574.to_numpy();y=d.GSE245386.to_numpy()
        same=(np.sign(x)==np.sign(y))&(x!=0)&(y!=0)
        ax.scatter(x,y,c=np.where(same,BLUE,ORANGE),s=9,alpha=.65,linewidths=0)
        ax.axhline(0,color=GRAY,lw=.6);ax.axvline(0,color=GRAY,lw=.6)
        rho=summary[(summary.config=='primary')&(summary.metric==metric)&(summary.statistic=='spearman_rho')].observed.iloc[0]
        ax.set_title(('Expression coavailability' if k==0 else 'LIANA 1.10.0 native priority')+f'\nSame 187 candidates; rho = {rho:.3f}')
        ax.set_xlabel('GSE174574 MCAO - Sham effect');ax.set_ylabel('GSE245386 MCAO - Sham effect')
        label(ax,'AB'[k])
    ax=axes[1,0]
    shown=['custom_score','lr_means','expr_prod','lrscore','magnitude_priority']
    labels=['Coavailability','LR mean','Expression product','SingleCellSignalR','Consensus priority']
    for ci,config in enumerate(['primary','reference_singlet']):
        for mi,metric in enumerate(shown):
            r=summary[(summary.config==config)&(summary.metric==metric)&(summary.statistic=='same_direction_all_fraction')].iloc[0]
            y=mi+(ci-.5)*.24
            color=BLUE if ci==0 else ORANGE
            ax.plot([r.null_q025*100,r.null_q975*100],[y,y],color=color,alpha=.42,lw=2.5)
            ax.scatter(r.null_median*100,y,marker='|',c=color,s=35)
            ax.scatter(r.observed*100,y,marker='o' if ci==0 else '^',c=color,s=24,zorder=3)
    ax.set_yticks(range(5),labels);ax.invert_yaxis();ax.set_xlim(0,100)
    ax.set_xlabel('Same direction / all common candidates (%)')
    ax.set_title('Observed points and 95% null reference intervals')
    ax.axhline(3.5,color=GRAY,lw=.6,ls=':')
    handles=[Line2D([],[],color=BLUE,marker='o',lw=0,label='Primary; 187'),Line2D([],[],color=ORANGE,marker='^',lw=0,label='Reference singlets; 174')]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,-.18),frameon=False);label(ax,'C')
    ax=axes[1,1]
    d=dyn[dyn.config=='primary']
    r=ds[(ds.config=='primary')&(ds.statistic=='same_direction_all_fraction')].iloc[0]
    assert len(d)==200
    ax.hist(d.same_direction_all_fraction*100,bins=np.linspace(0,100,21),color='#BBBBBB',edgecolor='white')
    ax.axvline(r.observed*100,color=BLUE,lw=1.5,label=f'Observed {r.observed*100:.1f}%')
    ax.axvline(r.null_median*100,color='black',ls='--',lw=.8,label='Null median 50%')
    ax.set_xlim(0,100);ax.set_xlabel('Concordance with gate recomputed (%)')
    ax.set_ylabel('Joint animal-label allocations')
    ax.set_title('Original 10% gate; exhaustive 20 x 10 null\nExact P = .005; six-test Holm P = .030')
    ax.legend(loc='upper left',bbox_to_anchor=(0,-.18),frameon=False);label(ax,'D')
    effects.to_csv(DATA/'figure2_common_effects.tsv.gz',sep='\t',index=False)
    summary.to_csv(DATA/'figure2_common_null_summary.tsv',sep='\t',index=False)
    d.to_csv(DATA/'figure2_original_gate_null.tsv',sep='\t',index=False)
    save(fig,'Figure2_expression_priority_and_null')

if __name__=='__main__':
    figure1();figure2()
