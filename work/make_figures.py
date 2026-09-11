from pathlib import Path
import json,string
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/figures';OUT.mkdir(parents=True,exist_ok=True)
DATA=ROOT.parent/'outputs/figure_source_data';DATA.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':400})
BLUE='#0072B2';ORANGE='#D55E00';GRAY='#8A8A8A';GREEN='#009E73'
names={'GSE174574':'Discovery (3 sham / 3 MCAO)','GSE245386':'Validation (3 sham / 2 MCAO)','GSE163752':'Sorted (6 paired pools)'}
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight')
 plt.close(fig);print(name,flush=True)
def label(ax,s):ax.text(-.12,1.06,s,transform=ax.transAxes,fontweight='bold',fontsize=10)
def points(ax,frame,value='score',xpos=(0,1),size=24):
 for j,cond in enumerate(['Sham','MCAO']):
  y=frame.loc[frame.condition==cond,value].to_numpy();off=np.linspace(-.11,.11,len(y))
  ax.scatter(np.full(len(y),xpos[j])+off,y,c=BLUE if j==0 else ORANGE,marker='o' if j==0 else '^',s=size,edgecolor='white',linewidth=.35,zorder=3)
  ax.plot([xpos[j]-.17,xpos[j]+.17],[y.mean()]*2,color='black',lw=1)
 return ax
# Figure 1: individual sample counts and lineage marker diagnosis.
fig,axes=plt.subplots(2,2,figsize=(7.2,5.5),layout='constrained',gridspec_kw={'height_ratios':[1,1.15]})
countall=[]
for k,acc in enumerate(names.keys()):
 if acc=='GSE163752':continue
 p=ROOT/'processed'/acc;c=pd.read_csv(p/'sample_cell_counts.tsv',sep='\t');c=c[(c.config=='primary')&c.cell_type.isin(['Astrocyte','Endothelial'])];countall.append(c)
 ax=axes[0,k]
 for j,typ in enumerate(['Astrocyte','Endothelial']):points(ax,c[c.cell_type==typ],value='n_cells',xpos=(j*3,j*3+1))
 ax.set_xticks([0,1,3,4],['Sham','MCAO','Sham','MCAO']);ax.set_ylim(bottom=0);ax.set_ylabel('Retained cells per biological sample');ax.set_title(names[acc]);ax.text(.18,-.24,'Astrocytes',transform=ax.transAxes);ax.text(.65,-.24,'Endothelium',transform=ax.transAxes);label(ax,string.ascii_uppercase[k])
 markers=['Slc1a2','Slc1a3','Aldh1l1','Pecam1','Cdh5','Kdr','Ptprc','Tyrobp','C1qa','Lyz2']
 d=pd.read_csv(ROOT/'results/identity_marker_audit.tsv',sep='\t');d=d[(d.dataset==acc)&(d.config=='primary')];t=d.groupby(['cell_type','gene']).fraction.mean().unstack('gene').reindex(index=['Astrocyte','Endothelial'],columns=markers)
 ax=axes[1,k];im=ax.imshow(t.to_numpy(),vmin=0,vmax=1,cmap='cividis',aspect='auto');ax.set_xticks(range(len(markers)),markers,rotation=60,ha='right');ax.set_yticks([0,1],['Astrocyte','Endothelial']);ax.set_title('Marker detection averaged across samples');label(ax,string.ascii_uppercase[k+2]);fig.colorbar(im,ax=ax,shrink=.65,label='Fraction of cells with >0 UMI')
 t.to_csv(DATA/f'figure1_markers_{acc}.tsv',sep='\t')
pd.concat(countall).to_csv(DATA/'figure1_sample_counts.tsv',sep='\t',index=False);save(fig,'Figure1_sample_and_identity_audit')
# Figure 2: all eligible effects, no cherry-picked scatter cloud.
fig,axes=plt.subplots(1,3,figsize=(7.2,2.75),layout='constrained')
summary=json.loads((ROOT/'results/cross_cohort/summary.json').read_text())
for k,typ in enumerate(['Astrocyte','Endothelial']):
 d=pd.read_csv(ROOT/f'results/cross_cohort/{typ}__gene_replication.tsv.gz',sep='\t');ax=axes[k]
 ax.scatter(d.log2FoldChange_discovery,d.log2FoldChange_validation,s=2,c=GRAY,alpha=.25,linewidths=0)
 sel=d[d.both_selected&d.same_direction];ax.scatter(sel.log2FoldChange_discovery,sel.log2FoldChange_validation,s=3,c=BLUE,alpha=.5,linewidths=0)
 ax.axhline(0,c='black',lw=.5);ax.axvline(0,c='black',lw=.5);ax.set_title(f'{typ}\nSpearman rho = {summary[typ]["rho_all"]:.2f}');ax.set_xlabel('Discovery log2 fold change');ax.set_ylabel('Validation log2 fold change');label(ax,string.ascii_uppercase[k]);d.to_csv(DATA/f'figure2_{typ}.tsv.gz',sep='\t',index=False)
d=pd.read_csv(ROOT/'results/cross_cohort/communication_replication.tsv.gz',sep='\t');d=d[d.both_eligible10];ax=axes[2]
ax.scatter(d.score_difference_discovery,d.score_difference_validation,s=8,c=np.where(d.same_direction,BLUE,ORANGE),alpha=.7,linewidths=0);ax.axhline(0,c='black',lw=.5);ax.axvline(0,c='black',lw=.5);ax.set_xlabel('Discovery score difference');ax.set_ylabel('Validation score difference');ax.set_title(f'Candidate coavailability\nrho = {summary["communication"]["rho_all_common"]:.2f}; n = {len(d)}');label(ax,'C');d.to_csv(DATA/'figure2_communication.tsv',sep='\t',index=False);save(fig,'Figure2_cross_cohort_concordance')
# Figure 3: sample-level scores for biologically interpretable examples.
pairs=json.loads((ROOT/'results/integrated/illustrative_pairs.json').read_text());fig,axes=plt.subplots(2,2,figsize=(7.2,5.4),layout='constrained');allpoints=[]
for k,(sender,receiver,ligand,receptor) in enumerate(pairs):
 ax=axes.flat[k]
 for j,acc in enumerate(['GSE174574','GSE245386']):
  d=pd.read_csv(ROOT/f'results/{acc}/primary__communication_sample_scores.tsv.gz',sep='\t');d=d[(d.sender==sender)&(d.receiver==receiver)&(d.ligand==ligand)&(d.receptor==receptor)].copy();d['dataset']=acc;allpoints.append(d);points(ax,d,xpos=(j*3,j*3+1),size=31)
 ax.set_xticks([0,1,3,4],['Sham\nn=3','MCAO\nn=3','Sham\nn=3','MCAO\nn=2']);ax.set_xlim(-.5,4.5);ax.set_ylim(bottom=0);ax.set_ylabel('Candidate coavailability score');direction='Astrocyte to endothelium' if sender=='Astrocyte' else 'Endothelium to astrocyte';ax.set_title(f'{ligand} - {receptor.replace("_", "/")}\n{direction}');ax.text(.08,-.3,'Discovery',transform=ax.transAxes);ax.text(.7,-.3,'Validation',transform=ax.transAxes);label(ax,string.ascii_uppercase[k])
pd.concat(allpoints).to_csv(DATA/'figure3_sample_scores.tsv',sep='\t',index=False);save(fig,'Figure3_candidate_sample_scores')
# Figure 4: gene-level barrier heterogeneity and sample panel summaries.
bbb=pd.read_csv(ROOT/'results/integrated/barrier_gene_effects.tsv',sep='\t');panels=json.loads((ROOT/'gene_panels.json').read_text());geneorder=sum([v['genes'] for v in panels.values()],[]);cohorts=['GSE174574','GSE245386','GSE163752']
fig=plt.figure(figsize=(7.2,6.6),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1.1,1]);ax=fig.add_subplot(gs[:,0]);mat=bbb.pivot(index='gene',columns='dataset',values='effect').reindex(index=geneorder,columns=cohorts);im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-3,vmax=3);ax.set_xticks([0,1,2],['Discovery','Validation','Sorted'],rotation=25,ha='right');ax.set_yticks(range(len(geneorder)),geneorder);ax.set_title('Endothelial transcript effects');label(ax,'A')
for i,g in enumerate(geneorder):
 for j,c in enumerate(cohorts):
  r=bbb[(bbb.gene==g)&(bbb.dataset==c)].iloc[0]
  if r.padj<.05:ax.scatter(j,i,marker='o',s=10,facecolors='none',edgecolors='black',linewidth=.7)
for y in [3.5,9.5,15.5]:ax.axhline(y,c='white',lw=1.2)
fig.colorbar(im,ax=ax,shrink=.45,label='Log2 effect (color clipped at +/-3)')
sc=pd.read_csv(ROOT/'results/integrated/endothelial_program_samples.tsv',sep='\t')
for k,panel in enumerate(['BBB_junction','BBB_transport']):
 ax=fig.add_subplot(gs[k,1]);sub=sc[sc.panel==panel]
 for j,acc in enumerate(cohorts[:2]):points(ax,sub[sub.dataset==acc],xpos=(3*j,3*j+1))
 ax.set_xticks([0,1,3,4],['Sham','MCAO','Sham','MCAO'],rotation=30,ha='right');ax.set_title('Junction-associated panel' if k==0 else 'Transport-associated panel');ax.set_ylabel('Mean log2 pseudobulk CPM+1');ax.text(.08,-.26,'Discovery',transform=ax.transAxes);ax.text(.69,-.26,'Validation',transform=ax.transAxes);label(ax,string.ascii_uppercase[k+1])
bbb.to_csv(DATA/'figure4_gene_effects.tsv',sep='\t',index=False);sc.to_csv(DATA/'figure4_program_scores.tsv',sep='\t',index=False);save(fig,'Figure4_barrier_transcription_heterogeneity')
# Figure 5: authentic coordinates, expression maps only; no invented anatomy.
maps={s:pd.read_csv(ROOT/f'results/spatial/{s}_spots.tsv.gz',sep='\t') for s in ['GSM7437221','GSM7437222']};genes=['Spp1','Itga5','Ptn','Ptprz1'];fig,axes=plt.subplots(2,4,figsize=(7.2,4.6),layout='constrained')
for j,g in enumerate(genes):
 vmax=max(np.percentile(np.concatenate([d.loc[d.qc_pass,g].to_numpy() for d in maps.values()]),99),.1)
 for i,(s,d) in enumerate(maps.items()):
  d=d[d.qc_pass];ax=axes[i,j];im=ax.scatter(d.imageX,d.imageY,c=d[g],s=2,cmap='viridis',vmin=0,vmax=vmax,linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([]);ax.spines[['bottom','left']].set_visible(False);ax.set_title(g)
  if j==0:ax.set_ylabel(('Control' if i==0 else 'Day 1')+f'\n{len(d):,} spots')
  if i==0:label(ax,string.ascii_uppercase[j])
 fig.colorbar(im,ax=axes[:,j],shrink=.45,orientation='horizontal',label='log1p(CP10k)',pad=.03)
for s,d in maps.items():d[['sample','time','barcode','qc_pass','imageX','imageY']+genes].to_csv(DATA/f'figure5_{s}.tsv.gz',sep='\t',index=False)
save(fig,'Figure5_spatial_expression_context')
# Supplement: held-out sorted constituent effects with unadjusted 95% confidence intervals.
c=pd.read_csv(ROOT/'results/integrated/illustrative_pair_components.tsv',sep='\t');fig,axes=plt.subplots(2,2,figsize=(7.2,5.6),layout='constrained')
for k,(sender,receiver,ligand,receptor) in enumerate(pairs):
 ax=axes.flat[k];d=c[(c.sender==sender)&(c.receiver==receiver)&(c.ligand==ligand)&(c.receptor==receptor)];gs=d[['role','gene']].drop_duplicates();labels=[]
 for i,r in enumerate(gs.itertuples(index=False)):
  labels.append(r.gene+' ('+r.role+')')
  for j,acc in enumerate(cohorts):
   z=d[(d.gene==r.gene)&(d.role==r.role)&(d.dataset==acc)]
   if not len(z):continue
   z=z.iloc[0];ax.errorbar(z.effect,i+(j-1)*.2,xerr=[[z.effect-z.ci_low],[z.ci_high-z.effect]],fmt=['o','^','s'][j],ms=3,color=[BLUE,ORANGE,GREEN][j],lw=.8,capsize=1)
 ax.axvline(0,color=GRAY,lw=.7);ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.set_title(ligand+' - '+receptor.replace('_','/'));ax.set_xlabel('Log2 expression effect (95% CI)');label(ax,string.ascii_uppercase[k])
handles=[Line2D([0],[0],marker=m,color=color,lw=.7,markersize=4,label=txt) for m,color,txt in zip(['o','^','s'],[BLUE,ORANGE,GREEN],['Discovery','Validation','Sorted; ipsi vs contra'])];fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False);c.to_csv(DATA/'figureS1_constituent_effects.tsv',sep='\t',index=False);save(fig,'FigureS1_constituent_expression')
# Supplement: all spatial sections, explicitly descriptive correlations.
d=pd.read_csv(ROOT/'results/spatial/descriptive_coexpression.tsv',sep='\t');d['pair']=d.ligand+' - '+d.receptor;order=['Control','Day 1','Day 3','Day 7a','Day 7b'];fig,axes=plt.subplots(1,2,figsize=(7.2,3),layout='constrained')
for k,col in enumerate(['spearman_same_spot','spearman_neighbor']):
 t=d.pivot(index='pair',columns='time',values=col).reindex(columns=order);ax=axes[k];im=ax.imshow(t,aspect='auto',vmin=-.6,vmax=.6,cmap='RdBu_r');ax.set_xticks(range(5),order,rotation=25,ha='right');ax.set_yticks(range(len(t)),t.index);ax.set_title('Same-spot expression' if k==0 else 'Six nearest spatial neighbors');label(ax,string.ascii_uppercase[k])
 for i in range(len(t)):
  for j in range(5):ax.text(j,i,f'{t.iloc[i,j]:.2f}',ha='center',va='center',fontsize=7)
 fig.colorbar(im,ax=ax,shrink=.6,label='Spearman rho; descriptive only')
d.to_csv(DATA/'figureS2_spatial_correlations.tsv',sep='\t',index=False);save(fig,'FigureS2_spatial_descriptive_correlations')
