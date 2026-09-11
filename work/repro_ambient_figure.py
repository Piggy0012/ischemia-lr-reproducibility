"""Figure 4: paired marker estimates, matched ranks, and gate coverage."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2';T=OUT/'tables';F=OUT/'figures';S=OUT/'figure_source_data'
for p in [F,S]:p.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':400})
BLUE='#0072B2';ORANGE='#D55E00';GRAY='#999999'

def label(ax,s):ax.text(-.13,1.06,s,transform=ax.transAxes,fontweight='bold',fontsize=11)

def main():
    m=pd.read_csv(T/'ambient_markers_per_animal.tsv',sep='\t')
    ranks=pd.read_csv(T/'ambient_rank_matched_summary.tsv',sep='\t')
    gate=pd.read_csv(T/'ambient_concordance_summary.tsv',sep='\t')
    dyn=pd.read_csv(T/'null_dynamic_gate_summary.tsv',sep='\t')
    assert m[['dataset','sample']].drop_duplicates().shape[0]==11
    fig,axes=plt.subplots(3,2,figsize=(7.2,7.5),layout='constrained',gridspec_kw={'height_ratios':[1.15,1.2,.85]})
    marker_rows=[];gate_rows=[]
    for j,typ in enumerate(['Astrocyte','Endothelial']):
        ax=axes[0,j]
        for k,gene in enumerate(['C1qa','Lyz2','Tyrobp']):
            d=m[(m.config=='primary')&(m.cell_type==typ)&(m.gene==gene)].sort_values(['dataset','sample'])
            assert len(d)==11;marker_rows.append(d)
            for i,r in enumerate(d.itertuples(index=False)):
                jitter=(i-5)*.016
                ys=np.log1p([r.raw_mean_count,r.corrected_mean_estimated_count])
                ax.plot([k-.13+jitter,k+.13+jitter],ys,color='#B8B8B8',lw=.5,zorder=1)
                ax.scatter(k-.13+jitter,ys[0],color=BLUE,s=10,zorder=2,linewidths=0)
                ax.scatter(k+.13+jitter,ys[1],color=ORANGE,s=11,marker='^',zorder=2,linewidths=0)
        ax.set_xticks(range(3),['C1qa','Lyz2','Tyrobp']);ax.set_ylim(bottom=0)
        ax.set_ylabel('log1p(mean counts per cell)');ax.set_title(typ+'; 11 paired animal libraries')
        label(ax,'AB'[j])
    metrics=['custom_score','lr_means','expr_prod','lrscore','native_magnitude_priority','diagnostic_unique_column_priority']
    labels=['Custom','LR mean','Product','LR score','Native','Unique']
    for j,config in enumerate(['primary','reference_singlet']):
        ax=axes[1,j]
        d=ranks[(ranks.config==config)&(ranks.comparison=='between_cohorts')]
        n=d.n_candidates.unique();assert len(n)==1
        for k,metric in enumerate(metrics):
            a=d[(d.metric==metric)&(d.count_source=='raw')].iloc[0]
            b=d[(d.metric==metric)&(d.count_source=='decontx')].iloc[0]
            assert a.n_candidates==b.n_candidates
            ax.plot([k-.12,k+.12],[a.spearman_rho,b.spearman_rho],color=GRAY,lw=.8,zorder=1)
            ax.scatter(k-.12,a.spearman_rho,color=BLUE,s=22,zorder=2)
            ax.scatter(k+.12,b.spearman_rho,color=ORANGE,s=26,marker='^',zorder=2)
        ax.axhline(0,color='#BBBBBB',lw=.6);ax.set_ylim(-1,1)
        ax.set_xticks(range(6),labels,rotation=35,ha='right');ax.set_ylabel('Cross-cohort effect Spearman rho')
        ax.set_title(('Primary' if j==0 else 'Reference-supported')+f'; same {n[0]} targets')
        label(ax,'CD'[j])
        ax=axes[2,j]
        raw=dyn[(dyn.config==config)&(dyn.statistic=='same_direction_all_fraction')]
        # The candidate count belongs to the observed gate, not the number
        # of allocations or a group of independent candidate observations.
        assert len(raw)==1
        nraw=int(raw.observed_n_candidates.iloc[0])
        counts=[nraw]
        for detection in ['gt0','ge1']:
            r=gate[(gate.config==config)&(gate.comparison=='between_cohort_corrected')&
                   (gate.corrected_detection==detection)&(gate.universe=='common_corrected_eligible10')]
            assert len(r)==1;counts.append(int(r.n_candidates.iloc[0]))
        descriptions=['Raw >0','Corrected >0','Corrected >=1']
        for k,(description,count) in enumerate(zip(descriptions,counts)):
            ax.barh(k,count,height=.55,color=[BLUE,ORANGE,'#999999'][k])
            ax.text(count+max(counts)*.025,k,str(count),va='center',fontsize=8)
            gate_rows.append({'config':config,'gate':description,'n_common_eligible':count})
        ax.set_yticks(range(3),descriptions);ax.invert_yaxis();ax.set_xlim(0,max(counts)*1.19)
        ax.set_xlabel('Common candidates passing 10% group gate')
        ax.set_title('Gate coverage; candidate sets may differ');label(ax,'EF'[j])
    fig.legend(handles=[Line2D([],[],marker='o',color=BLUE,linestyle='',label='Raw observed counts'),
        Line2D([],[],marker='^',color=ORANGE,linestyle='',label='DecontX estimated counts')],
        loc='outside lower center',ncol=2,frameon=False)
    for name,df in [('figure4_paired_markers.tsv',pd.concat(marker_rows)),('figure4_matched_rank_comparisons.tsv',ranks),('figure4_gate_coverage.tsv',pd.DataFrame(gate_rows))]:df.to_csv(S/name,sep='\t',index=False)
    fig.savefig(F/'Figure4_ambient_sensitivity.png',dpi=400)
    if '--export' in sys.argv:
        for ext in ['pdf','svg']:fig.savefig(F/f'Figure4_ambient_sensitivity.{ext}')
    plt.close(fig)
    print('Figure 4 preview generated'+('; vector export requested' if '--export' in sys.argv else ''))

if __name__=='__main__':main()
