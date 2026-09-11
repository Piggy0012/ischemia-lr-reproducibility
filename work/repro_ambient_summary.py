"""Descriptive ambient sensitivity; every gate and complete candidate set retained."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import stats
from repro_decontx_source import resolve
from revision_data import samples
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2/tables';OUT.mkdir(parents=True,exist_ok=True)
ACCS=['GSE174574','GSE245386'];CONFIGS=['primary','reference_singlet']
KEY=['sender','receiver','ligand','receptor'];TOL=1e-12
EXAMPLES=[('Astrocyte','Endothelial','Spp1','Itga5_Itgb1'),('Astrocyte','Endothelial','Timp3','Kdr'),
 ('Endothelial','Astrocyte','Ptn','Ptprz1'),('Endothelial','Astrocyte','Plat','Lrp1'),
 ('Astrocyte','Endothelial','Col4a1','Itga3_Itgb1')]

def compare(a,b):
    a=np.asarray(a,float).copy();b=np.asarray(b,float).copy()
    a[np.abs(a)<=TOL]=0;b[np.abs(b)<=TOL]=0
    nz=(a!=0)&(b!=0);same=nz&(np.sign(a)==np.sign(b))
    return {'n_candidates':len(a),'n_nonzero_both':int(nz.sum()),'n_same_direction':int(same.sum()),
        'n_zero_first':int((a==0).sum()),'n_zero_second':int((b==0).sum()),
        'same_direction_all_fraction':float(same.mean()) if len(a) else None,
        'same_direction_nonzero_fraction':float(same.sum()/nz.sum()) if nz.any() else None,
        'spearman_rho':float(stats.spearmanr(a,b).statistic) if len(a)>1 and np.std(a)>0 and np.std(b)>0 else None}

def main():
    coverage=[];comparisons=[];effects=[];examples=[];summary={}
    for config in CONFIGS:
        raw={};corrected={}
        for acc in ACCS:
            base=ROOT/('results' if config=='primary' else 'revision_analysis/results')/acc
            raw[acc]=pd.read_csv(base/f'{config}__communication.tsv.gz',sep='\t',float_precision='round_trip').set_index(KEY)
            corrected[acc]=pd.read_csv(ROOT/f'repro_decontx_analysis/results/{acc}/{config}__communication.tsv.gz',sep='\t',float_precision='round_trip').set_index(KEY)
            for key in EXAMPLES:
                a=raw[acc].loc[key] if key in raw[acc].index else None
                b=corrected[acc].loc[key] if key in corrected[acc].index else None
                examples.append(dict(zip(KEY,key))|{'dataset':acc,'config':config,
                    'raw_score_difference':None if a is None else a.score_difference,
                    'corrected_score_difference':None if b is None else b.score_difference,
                    'raw_eligible10':False if a is None else bool(a.eligible_10),
                    'corrected_eligible10_gt0':False if b is None else bool(b.eligible_10_gt0),
                    'corrected_eligible10_ge1':False if b is None else bool(b.eligible_10_ge1)})
        for gate in ['gt0','ge1']:
            eligible={};raweig={}
            for acc in ACCS:
                a=raw[acc];b=corrected[acc]
                raweig[acc]=a.index[a.eligible_10]
                eligible[acc]=b.index[b['eligible_10_'+gate]]
                common=raweig[acc].intersection(eligible[acc])
                coverage.append({'dataset':acc,'config':config,'corrected_detection':gate,
                    'n_raw_eligible10':len(raweig[acc]),'n_corrected_eligible10':len(eligible[acc]),
                    'n_common':len(common),'n_lost':len(raweig[acc].difference(eligible[acc])),
                    'n_gained':len(eligible[acc].difference(raweig[acc]))})
                comparisons.append({'comparison':'within_cohort_raw_vs_corrected','first':acc+'_raw','second':acc+'_corrected',
                    'config':config,'corrected_detection':gate,'universe':'common_raw_corrected_eligible10',
                    **compare(a.loc[common,'score_difference'],b.loc[common,'score_difference'])})
            common=eligible[ACCS[0]].intersection(eligible[ACCS[1]])
            comparisons.append({'comparison':'between_cohort_corrected','first':ACCS[0],'second':ACCS[1],
                'config':config,'corrected_detection':gate,'universe':'common_corrected_eligible10',
                **compare(corrected[ACCS[0]].loc[common,'score_difference'],corrected[ACCS[1]].loc[common,'score_difference'])})
            matched=common.intersection(raweig[ACCS[0]]).intersection(raweig[ACCS[1]])
            for source,tables in [('raw',raw),('corrected',corrected)]:
                comparisons.append({'comparison':'between_cohort_'+source,'first':ACCS[0],'second':ACCS[1],
                    'config':config,'corrected_detection':gate,'universe':'common_raw_corrected_both_cohorts_eligible10',
                    **compare(tables[ACCS[0]].loc[matched,'score_difference'],tables[ACCS[1]].loc[matched,'score_difference'])})
            for acc in ACCS:
                d=corrected[acc].loc[common].reset_index()
                d['corrected_detection']=gate;d['comparison_universe']='common_corrected_both_cohorts_eligible10'
                effects.append(d)
    for name,rows in [('coverage',coverage),('concordance_summary',comparisons),('fixed_examples',examples)]:
        pd.DataFrame(rows).to_csv(OUT/f'ambient_{name}.tsv',sep='\t',index=False)
    pd.concat(effects,ignore_index=True).to_csv(OUT/'ambient_cross_cohort_effects.tsv.gz',sep='\t',index=False)
    markers=pd.concat([pd.read_csv(f,sep='\t') for f in (ROOT/'repro_decontx_analysis/processed').glob('GSE*/*_markers.tsv')],ignore_index=True)
    assert markers[['dataset','sample']].drop_duplicates().shape[0]==11
    markers.to_csv(OUT/'ambient_markers_per_animal.tsv',sep='\t',index=False)
    model_samples=[];model_types=[];selected_sources=[]
    def model_summary(b):
        q=b.contamination.quantile([.25,.5,.75,.95])
        return {'n_cells':len(b),'mean_estimated_contamination':float(b.contamination.mean()),
            'median_estimated_contamination':float(q.loc[.5]),'q25_estimated_contamination':float(q.loc[.25]),
            'q75_estimated_contamination':float(q.loc[.75]),'q95_estimated_contamination':float(q.loc[.95]),
            'raw_total':float(b.raw_total.sum()),'corrected_total':float(b.corrected_total.sum()),
            'removed_total_fraction':float(1-b.corrected_total.sum()/b.raw_total.sum())}
    for acc in ACCS:
        for sample in samples(acc):
            folder,selection=resolve(acc,sample)
            a=json.loads((folder/'audit.json').read_text(encoding='utf-8'))
            b=pd.read_csv(folder/'contamination.tsv.gz',sep='\t')
            cells=pd.read_csv(folder/'cells.tsv.gz',sep='\t')
            assert np.array_equal(b.barcode,cells.barcode)
            b['original_rule']=cells.cell_type.to_numpy()
            base={'dataset':acc,'sample':sample}
            selected_sources.append(base|selection)
            model_samples.append(base|model_summary(b)|{'selected_max_iter':selection['selected_max_iter'],
                'selected_converged':selection['selected_converged'],'decontX_version':a['decontX'],
                'elapsed_seconds':a['elapsed_seconds']})
            for label_source in ['whole_brain_broad','original_rule']:
                for label,g in b.groupby(label_source):
                    model_types.append(base|{'label_source':label_source,'label':label}|model_summary(g))
    pd.DataFrame(model_samples).to_csv(OUT/'ambient_sample_contamination_summary.tsv',sep='\t',index=False)
    pd.DataFrame(model_types).to_csv(OUT/'ambient_cell_type_contamination_summary.tsv',sep='\t',index=False)
    pd.DataFrame(selected_sources).to_csv(OUT/'ambient_selected_sources.tsv',sep='\t',index=False)
    audit={'n_libraries':11,'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'no_new_pvalues':True,'fixed_identity_selections':CONFIGS,'corrected_detection_sensitivities':['>0','>=1 estimated count'],
        'normalization':'Each supplied matrix uses its own all-gene cell total',
        'source_provenance':'Fractional DecontX count estimates, not observed UMI',
        'interpretation':'Candidate loss under model-based correction is sensitivity, not proof that original expression was ambient.'}
    audit['convergence_source_selection']=selected_sources
    (OUT/'ambient_summary_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(pd.DataFrame(comparisons).to_string(index=False))

if __name__=='__main__':main()
