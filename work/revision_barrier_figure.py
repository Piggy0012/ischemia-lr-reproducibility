"""Descriptive BBB-gene sensitivity figure; no cell-level inference."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/revised'
d=pd.read_csv(OUT/'tables/barrier_genes_cell_selection_sensitivity.tsv',sep='\t')
panels=json.loads((ROOT/'gene_panels.json').read_text());genes=[g for p in panels.values() for g in p['genes']]
configs=['primary','singlet','reference_singlet','reference_only']
labels=['Original rules','Doublet filtered','Reference matched','Reference only']
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(7.2,7.1),layout='constrained',sharey=True)
for ax,acc,title,letter in zip(axes,['GSE174574','GSE245386'],['Discovery','Validation'],'AB'):
    sub=d[d.dataset==acc]
    matrix=sub.pivot(index='gene',columns='config',values='log2FoldChange').reindex(index=genes,columns=configs)
    fdr=sub.pivot(index='gene',columns='config',values='padj').reindex_like(matrix)
    im=ax.imshow(matrix,vmin=-3,vmax=3,cmap='RdBu_r',aspect='auto')
    ys,xs=np.where(fdr.to_numpy()<.05);ax.scatter(xs,ys,facecolors='none',edgecolors='black',s=12,lw=.65)
    for cut in [3.5,9.5,15.5]:ax.axhline(cut,color='white',lw=1.5)
    ax.set_xticks(range(4),labels,rotation=48,ha='right');ax.set_yticks(range(len(genes)),genes)
    ax.set_title(title,pad=8);ax.text(-.15,1.025,letter,transform=ax.transAxes,weight='bold',fontsize=13)
    ax.tick_params(length=0);ax.set_xlabel('Cell selection')
cb=fig.colorbar(im,ax=axes,shrink=.62,pad=.04);cb.set_label('MCAO vs Sham log2 fold change')
fig.get_layout_engine().set(rect=(0,.05,1,1))
fig.text(.5,.015,'Open circle: within-cell-type transcriptome BH FDR < 0.05',ha='center',fontsize=8)
(OUT/'figures').mkdir(parents=True,exist_ok=True)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figures/FigureS4_barrier_selection_sensitivity.{ext}',dpi=400,bbox_inches='tight')
plt.close(fig)
d.to_csv(OUT/'figure_source_data/figureS4_barrier_effects.tsv',sep='\t',index=False)
print('FigureS4_barrier_selection_sensitivity')
