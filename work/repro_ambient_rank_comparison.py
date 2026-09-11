"""Matched raw/corrected, all-animal candidate comparison; descriptive only."""
from pathlib import Path
import json, hashlib
import numpy as np,pandas as pd
from scipy import stats
from revision_data import samples
ROOT=Path(__file__).resolve().parent
TABLES=ROOT.parent/'outputs/reproducibility_v2/tables'
KEY=['sender','receiver','ligand','receptor']
RENAMES=dict(zip(['source','target','ligand_complex','receptor_complex'],KEY))
ACCS=['GSE174574','GSE245386'];CONFIGS=['primary','reference_singlet']
METRICS=['custom_score','lr_means','expr_prod','lrscore','native_magnitude_priority','diagnostic_unique_column_priority']
TOL=1e-12

def measure(x,y):
    x=np.asarray(x,float).copy();y=np.asarray(y,float).copy()
    x[np.abs(x)<=TOL]=0;y[np.abs(y)<=TOL]=0
    nz=(x!=0)&(y!=0);same=nz&(np.sign(x)==np.sign(y))
    return {'n_candidates':len(x),'n_nonzero_both':int(nz.sum()),'n_same_direction':int(same.sum()),
        'n_zero_first':int((x==0).sum()),'n_zero_second':int((y==0).sum()),
        'same_direction_all_fraction':float(same.mean()) if len(x) else None,
        'same_direction_nonzero_fraction':float(same.sum()/nz.sum()) if nz.any() else None,
        'spearman_rho':float(stats.spearmanr(x,y).statistic) if len(x)>1 and x.std()>0 and y.std()>0 else None}

def main():
    summary=[];effects=[];animal=[];keys=[];coverage=[]
    for config in CONFIGS:
        matrices={};manifests={};common=None
        for source in ['raw','decontx']:
            for acc in ACCS:
                target=pd.concat([pd.read_csv(ROOT/f'repro_liana_rank_diagnostic/{source}/{acc}/{config}__{s}__target.tsv.gz',sep='\t',float_precision='round_trip') for s in samples(acc)],ignore_index=True).rename(columns=RENAMES)
                folder=ROOT/('results' if config=='primary' else 'revision_analysis/results')/acc if source=='raw' else ROOT/'repro_decontx_analysis/results'/acc
                custom=pd.read_csv(folder/f'{config}__communication_sample_scores.tsv.gz',sep='\t',float_precision='round_trip').rename(columns={'score':'custom_score'})
                m=target[['sample','condition']].drop_duplicates().sort_values('sample')
                assert m['sample'].tolist()==sorted(samples(acc)) and m['sample'].is_unique
                manifests[acc]=m
                for metric in METRICS:
                    data=custom if metric=='custom_score' else target
                    mat=data.pivot(index=KEY,columns='sample',values=metric).reindex(columns=m['sample'])
                    complete=mat.index[np.isfinite(mat.to_numpy()).all(axis=1)]
                    common=complete if common is None else common.intersection(complete)
                    matrices[source,acc,metric]=mat
                    coverage.append({'config':config,'count_source':source,'dataset':acc,'metric':metric,'n_all_animal_complete':len(complete)})
        common=common.sort_values();assert len(common)>0
        k=common.to_frame(index=False);k.insert(0,'config',config);keys.append(k)
        effect={}
        for (source,acc,metric),mat in matrices.items():
            mat=mat.loc[common];case=manifests[acc].condition.eq('MCAO').to_numpy()
            x=mat.to_numpy();e=x[:,case].mean(axis=1)-x[:,~case].mean(axis=1)
            e[np.abs(e)<=TOL]=0;effect[source,acc,metric]=e
            d=common.to_frame(index=False);d['config']=config;d['count_source']=source;d['dataset']=acc;d['metric']=metric;d['mcao_minus_sham']=e;effects.append(d)
            d=mat.rename_axis(columns='sample').stack().rename('value').reset_index().merge(manifests[acc],on='sample',validate='many_to_one')
            d['config']=config;d['count_source']=source;d['dataset']=acc;d['metric']=metric;animal.append(d)
        for metric in METRICS:
            for source in ['raw','decontx']:
                summary.append({'config':config,'comparison':'between_cohorts','count_source':source,'metric':metric,
                    'first':ACCS[0],'second':ACCS[1],**measure(effect[source,ACCS[0],metric],effect[source,ACCS[1],metric])})
            for acc in ACCS:
                summary.append({'config':config,'comparison':'raw_vs_corrected_within_cohort','count_source':'paired','metric':metric,
                    'first':acc+'_raw','second':acc+'_decontx',**measure(effect['raw',acc,metric],effect['decontx',acc,metric])})
    for suffix,df in [('summary.tsv',pd.DataFrame(summary)),('effects.tsv.gz',pd.concat(effects)),
        ('animal_scores.tsv.gz',pd.concat(animal)),('candidate_keys.tsv',pd.concat(keys)),('coverage.tsv',pd.DataFrame(coverage))]:
        df.to_csv(TABLES/('ambient_rank_matched_'+suffix),sep='\t',index=False)
    (TABLES/'ambient_rank_matched_audit.json').write_text(json.dumps({'status':'complete','metrics':METRICS,
        'n_config_library_combinations':22,'raw_and_corrected_fixed_universe':True,'new_pvalues':False,
        'detection_rule':'Native LIANA >0 at expr_prop=0.1; >=1 custom gate analyzed separately',
        'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_selection':'default if converged, otherwise uniform maxIter=2000 numerical follow-up'},indent=2),encoding='utf-8')
    print(pd.DataFrame(summary).query("comparison=='between_cohorts'").to_string(index=False))

if __name__=='__main__':main()
