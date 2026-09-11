"""Aggregate DecontX estimates and compare expression candidates in fixed cells."""
from pathlib import Path
import os,sys,json,gc,hashlib
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import numpy as np
import pandas as pd
from scipy import sparse,stats
from revision_data import load_sample,samples
from repro_decontx_source import resolve

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'repro_decontx_analysis'
TABLES=ROOT.parent/'outputs/reproducibility_v2/tables'
CONFIGS=['primary','reference_singlet']
TYPES=['Astrocyte','Endothelial']
KEY=['sender','receiver','ligand','receptor']
MARKERS=['Slc1a2','Slc1a3','Aldh1l1','Pecam1','Cdh5','Kdr','Ptprc','Tyrobp','C1qa','Lyz2','Spp1','Timp3','Ptn','Plat','Ptprz1','Lrp1','Itga5','Itgb1','Col4a1','Itga3']

def process(acc,sample):
    folder=OUT/'processed'/acc;folder.mkdir(parents=True,exist_ok=True)
    done=folder/(sample+'_audit.json')
    source,selection_audit=resolve(acc,sample)
    if done.exists():
        prior=json.loads(done.read_text(encoding='utf-8'))
        if prior['source_audit_sha256']==selection_audit['selected_audit_sha256']:
            print('AMBIENT AGG EXISTS',sample,flush=True);return
    raw,q,genes=load_sample(acc,sample)
    corrected=sparse.load_npz(source/'corrected_counts.npz').tocsr()
    exported=pd.read_csv(source/'cells.tsv.gz',sep='\t')
    exported_genes=pd.read_csv(source/'genes.tsv.gz',sep='\t')['symbol']
    assert raw.shape==corrected.shape and np.array_equal(exported.barcode,q.barcode)
    assert np.array_equal(exported_genes.to_numpy(),genes.to_numpy())
    assert np.isfinite(corrected.data).all() and (corrected.data>=0).all()
    assert np.allclose(np.asarray(raw.sum(axis=1)).ravel(),q.n_umis)
    totals=np.asarray(corrected.sum(axis=1)).ravel()
    assert (totals>0).all(),f'{sample}: nonpositive corrected totals need explicit policy'
    identity=pd.read_csv(ROOT/f'revision_results/identity/{acc}/{sample}_identity.tsv.gz',sep='\t')
    assert np.array_equal(identity.barcode,q.barcode)
    condition='Sham' if ('sham' in q.prefix.iloc[0].lower() or '_WTC' in q.prefix.iloc[0]) else 'MCAO'
    tables={key:{'gene':genes} for key in ['estimated_pseudobulk','mean_logcp10k','fractions','fractions_ge1']}
    records=[];markers=[]
    for config in CONFIGS:
        for typ in TYPES:
            mask=(q.cell_type==typ).to_numpy()
            if config=='reference_singlet':
                mask&=(~identity.predicted_doublet.to_numpy(bool))&(identity.whole_brain_broad.to_numpy()==typ)&(identity.whole_brain_probability.to_numpy()>=.5)
            n=int(mask.sum());assert n>=30
            key=config+'__'+typ
            y=corrected[mask].astype(np.float64)
            tables['estimated_pseudobulk'][key]=np.asarray(y.sum(axis=0)).ravel()
            tables['fractions'][key]=np.asarray((y>0).mean(axis=0)).ravel()
            tables['fractions_ge1'][key]=np.asarray((y>=1).mean(axis=0)).ravel()
            yn=y.multiply((1e4/totals[mask])[:,None]).tocsr();yn.data=np.log1p(yn.data)
            tables['mean_logcp10k'][key]=np.asarray(yn.mean(axis=0)).ravel()
            records.append({'dataset':acc,'sample':sample,'condition':condition,'config':config,'cell_type':typ,'n_cells':n,'total_estimated_counts':float(y.sum())})
            raw_y=raw[mask].astype(np.float64)
            raw_fraction=np.asarray((raw_y>0).mean(axis=0)).ravel()
            raw_means=np.asarray(raw_y.mean(axis=0)).ravel()
            raw_y=raw_y.multiply((1e4/q.n_umis.to_numpy()[mask])[:,None]).tocsr();raw_y.data=np.log1p(raw_y.data)
            raw_log=np.asarray(raw_y.mean(axis=0)).ravel()
            for g in MARKERS:
                if g not in genes:continue
                gi=genes.get_loc(g)
                markers.append({'dataset':acc,'sample':sample,'condition':condition,'config':config,'cell_type':typ,'gene':g,
                    'raw_detection':raw_fraction[gi],'corrected_detection_gt0':tables['fractions'][key][gi],
                    'corrected_detection_ge1':tables['fractions_ge1'][key][gi],
                    'raw_mean_count':raw_means[gi],'corrected_mean_estimated_count':tables['estimated_pseudobulk'][key][gi]/n,
                    'raw_mean_logcp10k':raw_log[gi],'corrected_mean_logcp10k':tables['mean_logcp10k'][key][gi]})
    for suffix,d in tables.items():pd.DataFrame(d).to_csv(folder/f'{sample}_{suffix}.tsv.gz',sep='\t',index=False)
    pd.DataFrame(records).to_csv(folder/f'{sample}_sample_counts.tsv',sep='\t',index=False)
    pd.DataFrame(markers).to_csv(folder/f'{sample}_markers.tsv',sep='\t',index=False)
    done.write_text(json.dumps({'sample':sample,'dataset':acc,'shape':list(raw.shape),'all_positive_corrected_totals':True,
        'corrected_source_selection':selection_audit,
        'source_audit_sha256':hashlib.sha256((source/'audit.json').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'count_provenance':'Fractional DecontX estimates; no DESeq2 on these estimates',
        'normalization':'Corrected counts normalised by corrected total across all genes in the same cell',
        'fixed_cell_selections':CONFIGS},indent=2),encoding='utf-8')
    print('AMBIENT AGG DONE',sample,flush=True)
    del raw,corrected,y,yn,raw_y;gc.collect()

def summarize(acc):
    folder=OUT/'processed'/acc
    ms=samples(acc)
    count=pd.concat([pd.read_csv(folder/f'{s}_sample_counts.tsv',sep='\t') for s in ms],ignore_index=True)
    count.to_csv(folder/'sample_cell_counts.tsv',sep='\t',index=False)
    lr=pd.read_csv(ROOT/'literature/mouseconsensus.csv').drop_duplicates(['source_genesymbol','target_genesymbol'])
    for config in CONFIGS:
        case=count[count.config==config].drop_duplicates('sample').set_index('sample').loc[ms,'condition'].eq('MCAO').to_numpy()
        means={};fractions={};fractions1={}
        for typ in TYPES:
            col=config+'__'+typ
            for dest,suffix in [(means,'mean_logcp10k'),(fractions,'fractions'),(fractions1,'fractions_ge1')]:
                dest[typ]=pd.concat([pd.read_csv(folder/f'{s}_{suffix}.tsv.gz',sep='\t',index_col=0)[col].rename(s) for s in ms],axis=1)
                assert dest[typ].notna().all().all()
        rec=[];points=[]
        for sender,receiver in [TYPES,TYPES[::-1]]:
            for ligand,receptor in lr[['source_genesymbol','target_genesymbol']].itertuples(index=False,name=None):
                lg=ligand.split('_');rg=receptor.split('_')
                if not all(g in means[sender].index for g in lg+rg):continue
                score=np.sqrt(means[sender].loc[lg].min(axis=0).to_numpy()*means[receiver].loc[rg].min(axis=0).to_numpy())
                base=dict(zip(KEY,[sender,receiver,ligand,receptor]))|{'dataset':acc,'config':config}
                r=base|{'score_difference':float(score[case].mean()-score[~case].mean())}
                for source,label in [(fractions,'gt0'),(fractions1,'ge1')]:
                    ld=source[sender].loc[lg].min(axis=0).to_numpy();rd=source[receiver].loc[rg].min(axis=0).to_numpy()
                    r['eligible_10_'+label]=bool(max((ld[case]>=.1).sum(),(ld[~case]>=.1).sum())>=2 and max((rd[case]>=.1).sum(),(rd[~case]>=.1).sum())>=2)
                rec.append(r)
                for i,s in enumerate(ms):points.append(base|{'sample':s,'condition':'MCAO' if case[i] else 'Sham','score':score[i]})
        dest=OUT/'results'/acc;dest.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(rec).to_csv(dest/f'{config}__communication.tsv.gz',sep='\t',index=False)
        pd.DataFrame(points).to_csv(dest/f'{config}__communication_sample_scores.tsv.gz',sep='\t',index=False)
    print('AMBIENT CUSTOM DONE',acc,flush=True)

if __name__=='__main__':
    accs=[sys.argv[1]] if len(sys.argv)>1 else ['GSE174574','GSE245386']
    for acc in accs:
        chosen=[sys.argv[2]] if len(sys.argv)>2 else samples(acc)
        for s in chosen:process(acc,s)
        if len(sys.argv)<=2:summarize(acc)
