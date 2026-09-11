"""Read-only independent checks of controlled CellChat mapping and numerical comparators."""
from pathlib import Path
import collections
import gc
import gzip
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy import io, stats
from liana.method._pipe_utils._aggregate import _rank_aggregate

WORK=Path(__file__).resolve().parent
BASE=WORK/'repro_v3_cellchat'
PROBE=WORK/'repro_v3_cellchat_probe'
KEY=['source','target','ligand','receptor']
LM=['lr_means','expr_prod','lrscore']
sources={}


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def read(p,**kwargs):
    p=Path(p);sources[p.relative_to(WORK.parent).as_posix()]=sha(p)
    return pd.read_csv(p,sep='\t',float_precision='round_trip',**kwargs)


def main():
    native=read(PROBE/'CellChatDB.mouse_interaction.tsv',keep_default_na=False)
    complex_df=read(PROBE/'CellChatDB.mouse_complex.tsv',keep_default_na=False).set_index('resource_row_id')
    cofactors=read(PROBE/'CellChatDB.mouse_cofactor.tsv',keep_default_na=False).set_index('resource_row_id')
    def parts(value,tab):
        if not value:return []
        return [str(x) for x in tab.loc[value] if str(x)] if value in tab.index else [str(value)]
    def canonical(value,tab):return '_'.join(sorted(set(parts(value,tab))))
    native['cl']=native.ligand.map(lambda x:canonical(x,complex_df))
    native['cr']=native.receptor.map(lambda x:canonical(x,complex_df))
    protein=native[native.annotation.ne('Non-protein Signaling')].copy()
    multiplicity=protein.groupby(['cl','cr']).size()
    old=pd.read_csv(WORK/'literature/mouseconsensus.csv')
    oldcanon={('_'.join(sorted(set(a.split('_')))),'_'.join(sorted(set(b.split('_'))))) for a,b in zip(old.source_genesymbol,old.target_genesymbol)}
    manifest=read(WORK/'revision_sample_manifest.tsv')[['dataset','sample']]
    covered=None
    for acc,sample in manifest.itertuples(index=False,name=None):
        genes=set(read(WORK/f'revision_cache/{acc}/{sample}/genes.tsv.gz').symbol.astype(str))
        covered=genes if covered is None else covered & genes
    selected=[]
    for r in protein.itertuples(index=False):
        need=set(parts(r.ligand,complex_df)+parts(r.receptor,complex_df))
        for col in ['agonist','antagonist','co_A_receptor','co_I_receptor']:need.update(parts(getattr(r,col),cofactors))
        if (r.cl,r.cr) in oldcanon and multiplicity[r.cl,r.cr]==1 and need <= covered:selected.append(r.resource_row_id)
    controlled=read(BASE/'controlled_resource.tsv',keep_default_na=False)
    assert set(selected)==set(controlled.resource_row_id) and len(selected)==1548
    assert not controlled.duplicated(['canonical_ligand','canonical_receptor']).any()
    assert set(zip(controlled.canonical_ligand,controlled.canonical_receptor))==set(zip(protein.loc[protein.resource_row_id.isin(selected),'cl'],protein.loc[protein.resource_row_id.isin(selected),'cr']))
    fixed=read(BASE/'tables/fixed_complete_network_keys.tsv')
    target=read(BASE/'tables/fixed_complete_target_keys.tsv')
    common=None; full_counts=[]
    allowed=set(zip(controlled.canonical_ligand,controlled.canonical_receptor))
    for acc,sample in manifest.itertuples(index=False,name=None):
        p=WORK/f'repro_liana_rank_diagnostic/raw/{acc}/primary__{sample}__full_network.tsv.gz'
        d=read(p,usecols=['source','target','ligand_complex','receptor_complex']+LM)
        d['ligand']=d.ligand_complex.map(lambda x:'_'.join(sorted(set(x.split('_')))))
        d['receptor']=d.receptor_complex.map(lambda x:'_'.join(sorted(set(x.split('_')))))
        d=d[[(a,b) in allowed for a,b in zip(d.ligand,d.receptor)]]
        assert not d.duplicated(KEY).any() and np.isfinite(d[LM]).all().all()
        keys=set(map(tuple,d[KEY].to_numpy()));common=keys if common is None else common&keys
        full_counts.append({'sample':sample,'controlled_LIANA_network_rows':len(d)})
        del d;gc.collect()
    assert common==set(map(tuple,fixed.to_numpy())) and len(common)==1222
    target_keys={k for k in common if k[:2] in {('Astrocyte','Endothelial'),('Endothelial','Astrocyte')}}
    assert target_keys==set(map(tuple,target.to_numpy())) and len(target_keys)==32
    old187=read(WORK.parent/'outputs/reproducibility_v2/tables/null_fixed_candidate_keys.tsv')
    old187=old187[old187.config.eq('primary')].rename(columns={'sender':'source','receiver':'target'})
    oldkeys=set(map(tuple,old187[KEY].to_numpy()))
    nativepairs=set(zip(native.cl,native.cr))
    matched={k for k in oldkeys if k[2:] in nativepairs}
    assert len(oldkeys)==187 and matched==target_keys
    fullindex=pd.MultiIndex.from_frame(fixed);targetindex=pd.MultiIndex.from_frame(target)
    specs={'mean':('lr_means',False),'product':('expr_prod',False),'score':('lrscore',False)}
    dtype_rows=[];vectors={'float32':{},'float64':{}}
    conditions=read(BASE/'tables/sample_manifest.tsv')
    for sample in conditions['sample']:
        d=read(BASE/f'tables/{sample}__fixed_LIANA_components.tsv.gz').set_index(KEY).loc[fullindex]
        d32=d.copy();d32[LM]=d32[LM].astype(np.float32)
        a=1-_rank_aggregate(d32.reset_index(),specs,'rra')
        b=1-_rank_aggregate(d.reset_index(),specs,'rra')
        positions=fullindex.get_indexer(targetindex)
        vectors['float32'][sample]=a[positions];vectors['float64'][sample]=b[positions]
        dtype_rows.append({'sample':sample,'network_values_different':int((a!=b).sum()),
                           'max_absolute_network_difference':float(np.max(np.abs(a-b))),
                           'target_values_different':int((a[positions]!=b[positions]).sum())})
    dtype_stats=[]
    for dtype,values in vectors.items():
        effects=[]
        for acc in ['GSE174574','GSE245386']:
            m=conditions[conditions.dataset.eq(acc)];mat=np.column_stack([values[s] for s in m['sample']]);case=m.condition.eq('MCAO').to_numpy()
            e=mat[:,case].mean(axis=1)-mat[:,~case].mean(axis=1);e[np.abs(e)<=1e-12]=0;effects.append(e)
        a,b=effects;nz=(a!=0)&(b!=0);same=nz&(np.sign(a)==np.sign(b))
        dtype_stats.append({'input_dtype':dtype,'n_targets':32,'n_nonzero_both':int(nz.sum()),'n_same_direction':int(same.sum()),
                            'rho':float(stats.spearmanr(a,b).statistic)})
    # Directly verify 40 retained pilot cells using sparse memory-map row slices.
    pilot=BASE/'controlled/GSE245386/GSM7841720';cache=WORK/'revision_cache/GSE245386/GSM7841720'
    meta=read(pilot/'cells.tsv.gz');wanted=read(pilot/'genes.tsv').symbol.to_numpy()
    q=read(cache/'cells.tsv.gz');features=read(cache/'genes.tsv.gz').symbol.to_numpy()
    lookup={g:i for i,g in enumerate(wanted)};mapping=np.array([lookup.get(g,-1) for g in features])
    pos=pd.Index(q.barcode).get_indexer(meta.barcode);assert (pos>=0).all()
    dat=np.load(cache/'data.npy',mmap_mode='r');ind=np.load(cache/'indices.npy',mmap_mode='r');ptr=np.load(cache/'indptr.npy',mmap_mode='r')
    exported=io.mmread(pilot/'logcp10k_signal.mtx.gz').T.tocsr()
    cells=np.unique(np.linspace(0,len(meta)-1,40,dtype=int));maxerr=0
    for output_row in cells:
        rawrow=pos[output_row];sl=slice(ptr[rawrow],ptr[rawrow+1]);v=np.zeros(len(wanted),np.int32)
        indices=mapping[ind[sl]];keep=indices>=0
        np.add.at(v,indices[keep],dat[sl][keep])
        assert np.sum(dat[sl],dtype=np.int64)==meta.n_umis.iloc[output_row]==q.n_umis.iloc[rawrow]
        expected=np.log1p(v.astype(np.float32)*np.float32(1e4/meta.n_umis.iloc[output_row])).astype(np.float64)
        observed=exported.getrow(output_row).toarray().ravel();maxerr=max(maxerr,float(np.abs(expected-observed).max()))
        assert np.array_equal(expected,observed)
    assert np.array_equal(exported.data,exported.data.astype(np.float32).astype(np.float64))
    del exported,dat,ind,ptr;gc.collect()
    lossless=json.loads((pilot/'lossless_export_audit.json').read_text());telemetry=json.loads((pilot/'run_telemetry.json').read_text())
    assert lossless['source_rds_sha256']==sha(pilot/'native_network.rds')==telemetry['output_sha256']['native_network.rds']
    assert lossless['output_tsv_sha256']==sha(pilot/'full_network_lossless.tsv.gz')
    d=read(pilot/'full_network_lossless.tsv.gz')
    import rdata
    native_rds=rdata.read_rds(pilot/'native_network.rds')
    assert np.array_equal(np.asarray(native_rds['net']['prob']).ravel(order='F'),d.probability.to_numpy())
    assert np.array_equal(np.asarray(native_rds['net']['pval']).ravel(order='F'),d.pvalue.to_numpy())
    del native_rds;gc.collect()
    mapping=controlled.set_index('resource_row_id')[['canonical_ligand','canonical_receptor']]
    d['ligand']=d.interaction_name.map(mapping.canonical_ligand);d['receptor']=d.interaction_name.map(mapping.canonical_receptor)
    assert not d[KEY].isna().any().any() and not d.duplicated(KEY).any()
    cc=d.set_index(KEY).loc[fullindex];x=cc.probability.to_numpy();priority=(stats.rankdata(x,method='average')-1)/(len(x)-1)
    zero=x==0
    assert len(cc)==1222 and np.isfinite(x).all() and np.all(x>=0)
    nonzero=x>0
    assert not zero.any() or np.all(priority[zero]==(zero.sum()-1)/(2*(len(x)-1)))
    report={'status':'independent_checks_passed','original_187_exact_native_pairs_retained':len(matched),
            'original_187_not_exactly_in_native_resource':len(oldkeys-matched),
            'controlled_resource_rows':1548,'full_network_complete_keys':1222,'target_keys':32,
            'target_direction_counts':target.groupby(['source','target']).size().rename('n').reset_index().to_dict('records'),
            'target_unique_LR_pairs_ignoring_celltype_direction':len(target[['ligand','receptor']].drop_duplicates()),
            'n_targets_with_native_cofactor':sum(any(str(getattr(r,c))!='' for c in ['agonist','antagonist','co_A_receptor','co_I_receptor'])
                for k in target_keys for r in controlled[(controlled.canonical_ligand==k[2])&(controlled.canonical_receptor==k[3])].itertuples(index=False)),
            'controlled_resource_recomputed_independently':True,'fixed_full_network_recomputed_independently':True,
            'canonical_mapping_complete_for_exact_subunit_definition':True,'additional_target_loss_after_exact_resource_match':0,
            'per_library_eligible_rows':full_counts,'liana_dtype_sensitivity':dtype_rows,'liana_dtype_cross_cohort_sensitivity':dtype_stats,
            'pilot_normalization_cells_independently_checked':len(cells),'pilot_signal_genes_checked':len(wanted),
            'pilot_normalization_max_absolute_error':maxerr,'pilot_all_exported_positive_values_exact_float32_promotion':True,
            'pilot_normalization_check_scope':'40 sparse raw-cache cell rows across the retained cell order; every signaling gene; not every original library/cell',
            'pilot_probability_rows':len(d),'pilot_fixed_network_zero_strength':int(zero.sum()),
            'pilot_zero_strength_derived_percentile':float(priority[zero][0]) if zero.any() else None,
            'pilot_fixed_network_positive_strength_min':float(x[nonzero].min()) if nonzero.any() else None,
            'lossless_pilot_export_audit':lossless,'no_disease_CellChat_concordance_computed_by_review':True,
            'independent_Python_RDS_vs_TSV_probability_and_pvalue_bitwise_equal':True,
            'new_model_fit':False,'model_directory_files_modified':False,
            'review_conclusions':{
                'material_precision_issue':'Identified before CellChat cohort summary: restore saved LIANA component float32 before unique-column RRA. Analysis agent adopted the fix; observed 32-target concordance statistics unchanged in independent precision sensitivity.',
                'zero_rank_issue':'Average-rank percentiles have a library-dependent positive zero-strength floor. Final summary must retain direct strength and all-library-zero target diagnostics; analysis agent added them.',
                'final_model_completion_status':'Not certified by this review; native models and all-library lossless export completion are separate checks.',
                'inferential_scope':'Descriptive controlled-resource case study on 32 retained targets; no generic benchmark, independent biological replication or mechanism claim.',
                'final_provenance_suggestion':'Record fixed keys, manifests, LIANA component/custom means and lossless/model audits as summary input hashes, with summary output hashes.'},
            'source_sha256':sources,'review_script_sha256':sha(Path(__file__))}
    # Hash current code only after reading; it is not executed by this review.
    for name in ['repro_v3_independent_framework_plan.md','repro_v3_cellchat_run.py','repro_v3_cellchat_run.R','repro_v3_cellchat_summary.py','repro_v3_cellchat_export_lossless.R']:
        report['source_sha256']['work/'+name]=sha(WORK/name)
    (WORK/'repro_v3_cellchat_independent_review.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in {'source_sha256','per_library_eligible_rows','liana_dtype_sensitivity','lossless_pilot_export_audit'}},indent=2))


if __name__=='__main__':
    main()
