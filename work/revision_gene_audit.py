"""Marker audit of reference-added cells and BBB effects across cell selections."""
from pathlib import Path
import json,sys,gc
import numpy as np,pandas as pd
from revision_data import load_sample,samples
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/revised/tables';OUT.mkdir(parents=True,exist_ok=True)
MARKERS=['Aldh1l1','Slc1a3','Slc1a2','Gja1','Glul','Sox9','Aldoc','Pecam1','Cdh5','Esam','Tek','Erg','Vwf','Kdr','Ptprc','Tyrobp','Lyz2','C1qa','Rgs5','Pdgfrb','Acta2','Plp1','Mbp','Snap25']

def markers():
    records=[]
    for acc in ['GSE174574','GSE245386']:
        for sample in samples(acc):
            x,q,genes=load_sample(acc,sample)
            r=pd.read_csv(ROOT/f'revision_results/identity/{acc}/{sample}_identity.tsv.gz',sep='\t')
            assert np.array_equal(q.barcode,r.barcode)
            available=[g for g in MARKERS if g in genes]
            y=x[:,genes.get_indexer(available)].toarray();del x;gc.collect()
            norm=np.log1p(y*1e4/q.n_umis.to_numpy()[:,None])
            for typ in ['Astrocyte','Endothelial']:
                original=(q.cell_type==typ).to_numpy()
                reference=(r.whole_brain_broad==typ).to_numpy()&(r.whole_brain_probability.to_numpy()>=.5)&(~r.predicted_doublet.to_numpy(bool))
                for name,mask in [('original_reference_supported',original&reference),('reference_added',~original&reference),('original_not_retained',original&~reference)]:
                    for j,gene in enumerate(available):
                        records.append({'dataset':acc,'sample':sample,'cell_type':typ,'selection':name,'n_cells':int(mask.sum()),'gene':gene,'fraction_detected':float((y[mask,j]>0).mean()) if mask.any() else np.nan,'mean_logcp10k':float(norm[mask,j].mean()) if mask.any() else np.nan})
            print('MARKER AUDIT',sample,flush=True)
    pd.DataFrame(records).to_csv(OUT/'reference_added_cell_markers.tsv',sep='\t',index=False)

def bbb():
    panels=json.loads((ROOT/'gene_panels.json').read_text());records=[]
    for acc in ['GSE174574','GSE245386']:
        for config in ['primary','singlet','reference_singlet','reference_only']:
            parent=ROOT if config=='primary' else ROOT/'revision_analysis'
            d=pd.read_csv(parent/f'results/{acc}/{config}__Endothelial__de.tsv.gz',sep='\t').set_index('gene')
            for panel,desc in panels.items():
                for gene in desc['genes']:
                    row={'dataset':acc,'config':config,'panel':panel,'gene':gene,'tested':gene in d.index}
                    if gene in d.index:row.update(d.loc[gene].to_dict())
                    records.append(row)
    pd.DataFrame(records).to_csv(OUT/'barrier_genes_cell_selection_sensitivity.tsv',sep='\t',index=False)

if __name__=='__main__':
    if len(sys.argv)==1 or sys.argv[1]=='markers':markers()
    if len(sys.argv)==1 or sys.argv[1]=='bbb':bbb()
