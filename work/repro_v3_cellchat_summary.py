"""Descriptive same-universe comparison of native CellChat and cached LIANA scores."""
from pathlib import Path
import argparse, hashlib, json, gc
import numpy as np
import pandas as pd
from scipy import stats

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'repro_v3_cellchat';TABLES=OUT/'tables'
TABLES.mkdir(exist_ok=True)
KEY=['source','target','ligand','receptor'];TOL=1e-12
LMETRICS=['lr_means','expr_prod','lrscore']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(s):return '_'.join(sorted(set(str(s).split('_'))))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def read(p,**kw):return pd.read_csv(p,sep='\t',float_precision='round_trip',**kw)
def prepare():
    resource=read(OUT/'controlled_resource.tsv')
    pairs=set(zip(resource.canonical_ligand,resource.canonical_receptor))
    mapping=resource.set_index('resource_row_id')[['canonical_ligand','canonical_receptor','has_cofactor']]
    manifest=read(ROOT/'revision_sample_manifest.tsv')[['dataset','sample']]
    common=None;hashes={};counts=[];frames={}
    for acc,sample in manifest.itertuples(index=False,name=None):
        p=ROOT/f'repro_liana_rank_diagnostic/raw/{acc}/primary__{sample}__full_network.tsv.gz'
        d=read(p,usecols=['source','target','ligand_complex','receptor_complex','condition']+LMETRICS)
        d['ligand']=d.ligand_complex.map(canon);d['receptor']=d.receptor_complex.map(canon)
        d=d[[*KEY,'condition',*LMETRICS]]
        d=d[[(a,b) in pairs for a,b in zip(d.ligand,d.receptor)]]
        assert not d.duplicated(KEY).any() and d[LMETRICS].notna().all().all()
        ix=pd.MultiIndex.from_frame(d[KEY]);common=ix if common is None else common.intersection(ix)
        frames[sample]=d;hashes[str(p.relative_to(ROOT))]=sha(p)
        manifest.loc[manifest['sample']==sample,'condition']=d.condition.iloc[0]
        counts.append({'dataset':acc,'sample':sample,'shared_resource_LIANA_eligible_full_network_rows':len(d)})
    common=common.sort_values();assert len(common)>0
    target=common[((common.get_level_values('source')=='Astrocyte')&(common.get_level_values('target')=='Endothelial'))|
                  ((common.get_level_values('source')=='Endothelial')&(common.get_level_values('target')=='Astrocyte'))]
    assert len(target)>0
    common.to_frame(index=False).to_csv(TABLES/'fixed_complete_network_keys.tsv',sep='\t',index=False)
    target.to_frame(index=False).to_csv(TABLES/'fixed_complete_target_keys.tsv',sep='\t',index=False)
    manifest.to_csv(TABLES/'sample_manifest.tsv',sep='\t',index=False)
    pd.DataFrame(counts).to_csv(TABLES/'eligibility_attrition.tsv',sep='\t',index=False)
    for sample,d in frames.items():
        d.set_index(KEY).loc[common].to_csv(TABLES/(sample+'__fixed_LIANA_components.tsv.gz'),sep='\t')
    report={'status':'fixed_network_prepared','n_samples':len(manifest),'n_shared_resource_pairs':len(resource),
      'n_complete_full_network_keys':len(common),'n_complete_target_keys':len(target),
      'context_types_in_fixed_network':sorted(set(common.get_level_values('source'))|set(common.get_level_values('target'))),
      'condition_counts':manifest.groupby(['dataset','condition']).size().rename('n').reset_index().to_dict('records'),
      'source_LIANA_full_network_hashes':hashes,'resource_sha256':sha(OUT/'controlled_resource.tsv'),
      'plan_sha256':sha(ROOT/'repro_v3_independent_framework_plan.md'),
      'gate_note':'Complete score rows in all 11 LIANA libraries on the controlled resource, before looking at new CellChat scores. Conditions do not determine eligibility.',
      'no_biological_results_compared_yet':True}
    dump(TABLES/'fixed_universe_audit.json',report)
    print(json.dumps({k:report[k] for k in ['status','n_complete_full_network_keys','n_complete_target_keys','context_types_in_fixed_network']},indent=2))
    return common,target,manifest,mapping
def summarize():
    from liana.method._pipe_utils._aggregate import _rank_aggregate
    common=pd.MultiIndex.from_frame(read(TABLES/'fixed_complete_network_keys.tsv'))
    target=pd.MultiIndex.from_frame(read(TABLES/'fixed_complete_target_keys.tsv'))
    manifest=read(TABLES/'sample_manifest.tsv');resource=read(OUT/'controlled_resource.tsv')
    mapping=resource.set_index('resource_row_id')[['canonical_ligand','canonical_receptor','has_cofactor']]
    rra_specs={'CellPhoneDB':('lr_means',False),'Connectome':('expr_prod',False),'SingleCellSignalR':('lrscore',False)}
    rows=[];ties=[];input_hashes={};network_rows=[]
    scientific_input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [
        TABLES/'fixed_complete_network_keys.tsv',TABLES/'fixed_complete_target_keys.tsv',
        TABLES/'sample_manifest.tsv',TABLES/'fixed_universe_audit.json',OUT/'controlled_resource.tsv',
        OUT/'installed_package_verified.json',ROOT/'repro_v3_independent_framework_plan.md',
        ROOT/'repro_v3_cellchat_precision_amendment.md',ROOT/'repro_v3_cellchat_structural_audit_note.md']}
    for acc,sample,condition in manifest[['dataset','sample','condition']].itertuples(index=False,name=None):
        folder=OUT/f'controlled/{acc}/{sample}'
        a=json.loads((folder/'run_telemetry.json').read_text());assert a['status']=='complete'
        export_audit=json.loads((folder/'lossless_export_audit.json').read_text())
        assert export_audit['source_rds_sha256']==a['output_sha256']['native_network.rds']
        assert export_audit['exact_probability_roundtrip'] and export_audit['exact_pvalue_roundtrip']
        p=folder/'full_network_lossless.tsv.gz';assert export_audit['output_tsv_sha256']==sha(p)
        d=read(p);d['ligand']=d.interaction_name.map(mapping.canonical_ligand);d['receptor']=d.interaction_name.map(mapping.canonical_receptor)
        assert d[KEY].notna().all().all() and not d.duplicated(KEY).any()
        cc=d.set_index(KEY).loc[common];l=read(TABLES/(sample+'__fixed_LIANA_components.tsv.gz')).set_index(KEY).loc[common]
        # Archived component columns originated as float32. Restore their
        # scientific dtype before SciPy ranking/RRA, not merely decimal text.
        l[LMETRICS]=l[LMETRICS].astype(np.float32)
        assert np.isfinite(cc[['probability','pvalue']]).all().all()
        assert len(cc)==len(l)==len(common)
        cc_priority=(stats.rankdata(cc.probability.to_numpy(),method='average')-1)/(len(common)-1)
        liana_priority=1-_rank_aggregate(l.reset_index().copy(),rra_specs,'rra')
        full=pd.DataFrame({'cellchat_probability':cc.probability.to_numpy(),'cellchat_internal_pvalue':cc.pvalue.to_numpy(),
                           'cellchat_strength_percentile':cc_priority,'LIANA_unique_column_RRA_priority':liana_priority},index=common)
        for col in LMETRICS:full[col]=l[col].to_numpy()
        full.to_csv(TABLES/(sample+'__fixed_full_network_scores.tsv.gz'),sep='\t')
        chosen=full.loc[target].copy()
        mean=read(ROOT/f'processed/{acc}/{sample}_mean_logcp10k.tsv.gz',index_col=0)
        for source in [p,folder/'run_telemetry.json',folder/'model_audit.json',folder/'input_audit.json',
                       folder/'lossless_export_audit.json',TABLES/(sample+'__fixed_LIANA_components.tsv.gz'),
                       ROOT/f'processed/{acc}/{sample}_mean_logcp10k.tsv.gz']:
            scientific_input_hashes[str(source.relative_to(ROOT))]=sha(source)
        chosen['custom_coavailability']=[np.sqrt(mean.loc[ligand.split('_'),'primary__'+sender].min()*mean.loc[receptor.split('_'),'primary__'+receiver].min()) for sender,receiver,ligand,receptor in target]
        chosen['dataset']=acc;chosen['sample']=sample;chosen['condition']=condition
        rows.append(chosen.reset_index())
        ties.append({'dataset':acc,'sample':sample,'fixed_network_rows':len(common),'cellchat_zero_network_rows':int((cc.probability==0).sum()),
          'cellchat_unique_probability_values':cc.probability.nunique(),'cellchat_zero_target_rows':int((chosen.cellchat_probability==0).sum()),
          'cellchat_target_rows':len(target),'liana_unique_priority_values':int(pd.Series(liana_priority).nunique()),
          'cellchat_zero_strength_percentile_floor':float(cc_priority[cc.probability.to_numpy()==0][0]) if (cc.probability==0).any() else None,
          'probability_min_positive':float(cc.loc[cc.probability>0,'probability'].min()) if (cc.probability>0).any() else None})
        input_hashes[str(p.relative_to(ROOT))]=sha(p)
        del d,cc,l,mean,full;gc.collect()
    scores=pd.concat(rows,ignore_index=True)
    metrics=['custom_coavailability','cellchat_probability','cellchat_strength_percentile','LIANA_unique_column_RRA_priority']+LMETRICS
    effects=[];summaries=[]
    cofactor={(r.canonical_ligand,r.canonical_receptor):bool(r.has_cofactor) for _,r in resource.iterrows()}
    for metric in metrics:
        vectors=[]
        for acc in ['GSE174574','GSE245386']:
            m=manifest[manifest.dataset==acc];mat=scores[scores.dataset==acc].pivot(index=KEY,columns='sample',values=metric).loc[target,m['sample']]
            case=m.condition.to_numpy()=='MCAO';values=mat.to_numpy(dtype=float)
            raw=values[:,case].mean(axis=1)-values[:,~case].mean(axis=1)
            effect=raw.copy();effect[np.abs(effect)<=TOL]=0;vectors.append(effect)
            e=target.to_frame(index=False);e['dataset']=acc;e['metric']=metric;e['raw_mcao_minus_sham']=raw;e['mcao_minus_sham']=effect
            e['nonzero_effect_zeroed_at_1e_12']=(raw!=0)&(effect==0);effects.append(e)
        a,b=vectors;nz=(a!=0)&(b!=0);same=nz&(np.sign(a)==np.sign(b))
        scopes=[('all_fixed_targets',np.ones(len(target),bool)),
                ('Astrocyte_to_Endothelial',target.get_level_values('source')=='Astrocyte'),
                ('Endothelial_to_Astrocyte',target.get_level_values('source')=='Endothelial'),
                ('without_native_cofactor',np.array([not cofactor[(lig,rec)] for _,_,lig,rec in target])),
                ('with_native_cofactor',np.array([cofactor[(lig,rec)] for _,_,lig,rec in target]))]
        for scope,mask in scopes:
            n=int(mask.sum())
            rho=float(stats.spearmanr(a[mask],b[mask]).statistic) if n>1 and np.ptp(a[mask])>0 and np.ptp(b[mask])>0 else None
            summaries.append({'scope':scope,'metric':metric,'n_candidates':n,'n_nonzero_both':int((nz&mask).sum()),'n_same_direction':int((same&mask).sum()),
                'same_direction_all_fraction':float((same&mask).sum()/n) if n else None,'same_direction_nonzero_fraction':float((same&mask).sum()/(nz&mask).sum()) if (nz&mask).sum() else None,
                'n_zero_discovery':int(((a==0)&mask).sum()),'n_zero_external':int(((b==0)&mask).sum()),'spearman_rho':rho,'disease_inferential_pvalues_computed':False})
    scores.to_csv(TABLES/'target_sample_scores.tsv.gz',sep='\t',index=False)
    effect_table=pd.concat(effects,ignore_index=True)
    effect_table.to_csv(TABLES/'target_disease_effects.tsv.gz',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(TABLES/'cross_cohort_concordance.tsv',sep='\t',index=False)
    pd.DataFrame(ties).to_csv(TABLES/'tie_zero_audit.tsv',sep='\t',index=False)
    zeros=[]
    for key,g in scores.groupby(KEY,sort=True):
        z=dict(zip(KEY,key));z['zero_strength_library_count']=int((g.cellchat_probability==0).sum());z['all_11_strength_zero']=bool((g.cellchat_probability==0).all())
        for acc,label in [('GSE174574','discovery'),('GSE245386','external')]:
            sub=g[g.dataset==acc];case=sub.condition=='MCAO'
            z[label+'_zero_strength_library_count']=int((sub.cellchat_probability==0).sum())
            for col,labelmetric in [('cellchat_probability','direct_strength'),('cellchat_strength_percentile','derived_percentile')]:
                value=float(sub.loc[case,col].mean()-sub.loc[~case,col].mean())
                z[label+'_'+labelmetric+'_raw_effect']=value
                z[label+'_'+labelmetric+'_effect_at_1e_12']=0.0 if abs(value)<=TOL else value
        z['allzero_strength_with_nonzero_derived_disease_effect']=z['all_11_strength_zero'] and (z['discovery_derived_percentile_effect_at_1e_12']!=0 or z['external_derived_percentile_effect_at_1e_12']!=0)
        zeros.append(z)
    zero_table=pd.DataFrame(zeros);zero_table.to_csv(TABLES/'zero_strength_target_audit.tsv',sep='\t',index=False)
    audit={'status':'complete','arm':'controlled_shared_resource','n_network_keys':len(common),'n_target_keys':len(target),'n_samples':len(manifest),
      'input_CellChat_full_network_hashes':input_hashes,'script_sha256':sha(__file__),'frozen_plan_sha256':sha(ROOT/'repro_v3_independent_framework_plan.md'),
      'scientific_input_file_sha256':scientific_input_hashes,
      'rank_formula':'(ascending average rank of native CellChat probability - 1)/(fixed full-network N - 1); analyst-derived, not native CellChat consensus',
      'liana_unique_score_specs':rra_specs,'zero_effect_tolerance':TOL,'cofactor_subgroups_note':'Exploratory descriptive stratification; no independent tests or multiple-validation claim.',
      'liana_component_dtype_before_RRA':'float32, restored after round_trip TSV parsing; native CellChat probability remains float64',
      'disease_group_mean_arithmetic':'float64 for every metric',
      'all_11_zero_strength_targets':int(zero_table.all_11_strength_zero.sum()),
      'all_11_zero_strength_targets_with_nonzero_derived_disease_effect':int(zero_table.allzero_strength_with_nonzero_derived_disease_effect.sum()),
      'zero_rank_interpretation':'Average-rank percentiles assign the zero-strength tie block a positive floor that can differ by library. Nonzero derived rank effects of all-zero native strengths are relative-background behavior, not reproduced communication evidence.',
      'all_metrics_use_same_target_keys':True,'existing_LIANA_ranks_not_reused':True,'native_internal_pvalues_not_used_for_filtering_or_disease_inference':True,
      'new_inferential_pvalues_computed':False,'native_resource_arm_completed':False}
    exported=[TABLES/n for n in ['target_sample_scores.tsv.gz','target_disease_effects.tsv.gz','cross_cohort_concordance.tsv','tie_zero_audit.tsv','zero_strength_target_audit.tsv']]
    exported += [TABLES/(s+'__fixed_full_network_scores.tsv.gz') for s in manifest['sample']]
    audit['output_file_sha256']={str(p.relative_to(ROOT)):sha(p) for p in exported}
    dump(TABLES/'summary_audit.json',audit)
    print(pd.DataFrame(summaries).query("scope=='all_fixed_targets'").to_string(index=False))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prepare-only',action='store_true');args=parser.parse_args()
    if args.prepare_only:prepare()
    else:summarize()
