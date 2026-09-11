"""Original cohort pair restricted to the same 99 complete third-cohort edges."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parent
TABLES=ROOT.parent/'outputs/reproducibility_v2/tables'
PATH=TABLES/'third_rank_diagnostic_library_scores.tsv.gz'
KEY=['sender','receiver','ligand','receptor']
METRICS=['custom_score','lr_means','expr_prod','lrscore','native_magnitude_priority','diagnostic_unique_column_priority']
d=pd.read_csv(PATH,sep='\t',float_precision='round_trip')
d=d[(d.universe=='all_three_cohorts_complete') & (d.dataset!='GSE332910')]
assert not d.duplicated(['dataset','metric','sample']+KEY).any()
assert set(d.dataset)=={'GSE174574','GSE245386'}
rows=[]
all_effects=[]
for metric in METRICS:
    effects=[]
    common=None
    for acc,ncase in [('GSE174574',3),('GSE245386',2)]:
        b=d[(d.metric==metric)&(d.dataset==acc)]
        m=b[['sample','condition']].drop_duplicates().sort_values('sample')
        assert m.condition.value_counts().to_dict()=={'Sham':3,'MCAO':ncase}
        mat=b.pivot(index=KEY,columns='sample',values='value').reindex(columns=m['sample']).sort_index()
        assert len(mat)==99 and np.isfinite(mat.to_numpy()).all()
        assert common is None or common.equals(mat.index)
        common=mat.index
        case=m.condition.to_numpy()=='MCAO'
        x=mat.to_numpy()
        delta=x[:,case].mean(axis=1)-x[:,~case].mean(axis=1)
        delta[np.abs(delta)<=1e-12]=0
        effects.append(delta)
        e=common.to_frame(index=False)
        e['dataset'],e['metric'],e['mcao_minus_sham']=acc,metric,delta
        all_effects.append(e)
    a,b=effects
    nz=(a!=0)&(b!=0)
    same=nz&(np.sign(a)==np.sign(b))
    rows.append({'discovery_dataset':'GSE174574','validation_dataset':'GSE245386',
        'universe':'all_three_cohorts_complete','config':'primary','count_source':'raw','metric':metric,
        'n_candidates':99,'n_nonzero_both':int(nz.sum()),'n_same_direction':int(same.sum()),
        'same_direction_all_fraction':float(same.mean()),'same_direction_nonzero_fraction':float(same.sum()/nz.sum()),
        'n_zero_discovery':int((a==0).sum()),'n_zero_validation':int((b==0).sum()),
        'spearman_rho':float(stats.spearmanr(a,b).statistic),'inferential_pvalues_computed':False})
out=TABLES/'third_rank_diagnostic_original_pair_on99.tsv'
pd.DataFrame(rows).to_csv(out,sep='\t',index=False)
pd.concat(all_effects,ignore_index=True).to_csv(TABLES/'third_rank_diagnostic_original_pair_on99_effects.tsv.gz',sep='\t',index=False)
audit={'status':'complete','input_library_scores_sha256':hashlib.sha256(PATH.read_bytes()).hexdigest(),
    'same_preexisting_17_library_candidate_universe':True,'n_candidates':99,
    'old_summary_and_input_scores_not_modified':True,'tsv_float_parser':'round_trip',
    'no_new_inferential_pvalues':True,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'comparisons':rows}
(TABLES/'third_rank_diagnostic_original_pair_on99_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))
