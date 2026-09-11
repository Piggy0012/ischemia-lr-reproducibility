"""Uniform maxIter=2000 refits of default non-converged libraries; no candidate-based selection."""
from pathlib import Path
import os,json,sys,hashlib
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'repro_decontx'
EXT=ROOT/'repro_decontx_maxiter2000'
os.environ['DECONTX_OUTPUT_ROOT']=str(EXT)
os.environ['DECONTX_MAX_ITER']='2000'
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS']:
    os.environ[name]='1'
import numpy as np,pandas as pd
from scipy import sparse
from repro_decontx import run,summarize

def main():
    manifest=pd.read_csv(ROOT/'revision_sample_manifest.tsv',sep='\t')
    selected=[]
    for r in manifest.itertuples(index=False):
        p=BASE/r.dataset/r.sample/'audit.json'
        assert p.exists(),f'All 11 default outputs must finish first: {p}'
        a=json.loads(p.read_text())
        assert a['requested_parameters']['maxIter']==500
        if not a['convergence_threshold_reached']:selected.append((r.dataset,r.sample))
    EXT.mkdir(exist_ok=True)
    plan_bytes=(ROOT/'repro_decontx_convergence_amendment.md').read_bytes()
    driver_bytes=Path(__file__).read_bytes()
    (EXT/'analysis_plan_snapshot.md').write_bytes(plan_bytes)
    (EXT/'extension_driver_source_snapshot.py').write_bytes(driver_bytes)
    (EXT/'selection.json').write_text(json.dumps({'rule':'all default non-converged libraries; no candidate outcome selection',
        'maxIter':2000,'samples':selected,'default_samples':len(manifest),
        'amendment_sha256':hashlib.sha256(plan_bytes).hexdigest(),
        'extension_driver_sha256':hashlib.sha256(driver_bytes).hexdigest()},indent=2),encoding='utf-8')
    comparisons=[]
    for acc,sample in selected:
        run(acc,sample)
        base=BASE/acc/sample;ext=EXT/acc/sample
        a=json.loads((base/'audit.json').read_text());b=json.loads((ext/'audit.json').read_text())
        expected_parameters=dict(a['requested_parameters']);expected_parameters['maxIter']=2000
        assert b['requested_parameters']==expected_parameters
        assert a['R']==b['R'] and a['decontX']==b['decontX']
        assert a['input']['input_sha256']==b['input']['input_sha256']
        assert a['input']['identity_sha256']==b['input']['identity_sha256']
        assert np.array_equal(pd.read_csv(base/'genes.tsv.gz',sep='\t').symbol,pd.read_csv(ext/'genes.tsv.gz',sep='\t').symbol)
        c=pd.read_csv(base/'contamination.tsv.gz',sep='\t');d=pd.read_csv(ext/'contamination.tsv.gz',sep='\t')
        assert np.array_equal(c.barcode,d.barcode) and np.array_equal(c.whole_brain_broad,d.whole_brain_broad)
        x=sparse.load_npz(base/'corrected_counts.npz');y=sparse.load_npz(ext/'corrected_counts.npz')
        assert x.shape==y.shape
        if np.array_equal(x.indptr,y.indptr) and np.array_equal(x.indices,y.indices):
            diff=y.data-x.data
            l1=float(np.abs(diff).sum());maxdiff=float(np.max(np.abs(diff)))
        else:
            diff=(y-x).tocsr();l1=float(np.abs(diff.data).sum());maxdiff=float(np.max(np.abs(diff.data)))
        comparisons.append({'dataset':acc,'sample':sample,'default_last_iteration':a['last_logged_iteration'],
            'default_last_max_divergence':a['last_logged_max_divergence'],
            'extended_last_iteration':b['last_logged_iteration'],'extended_last_max_divergence':b['last_logged_max_divergence'],
            'extended_converged':b['convergence_threshold_reached'],'identical_inputs':True,
            'default_median_contamination':float(c.contamination.median()),'extended_median_contamination':float(d.contamination.median()),
            'mean_absolute_contamination_change':float(np.mean(np.abs(d.contamination-c.contamination))),
            'max_absolute_contamination_change':float(np.max(np.abs(d.contamination-c.contamination))),
            'count_matrix_l1_difference':l1,'count_matrix_l1_divided_by_raw_total':l1/float(c.raw_total.sum()),
            'largest_absolute_estimated_count_change':maxdiff})
        del x,y,diff
    summarize()
    pd.DataFrame(comparisons).to_csv(EXT/'convergence_comparison.tsv',sep='\t',index=False)
    (EXT/'completion.json').write_text(json.dumps({'selected_samples':len(selected),'completed_refits':len(comparisons),
        'all_completed':len(selected)==len(comparisons),'all_refits_converged':all(r['extended_converged'] for r in comparisons),
        'comparisons':comparisons},indent=2),encoding='utf-8')
    print(json.dumps(comparisons,indent=2),flush=True)

if __name__=='__main__':main()
