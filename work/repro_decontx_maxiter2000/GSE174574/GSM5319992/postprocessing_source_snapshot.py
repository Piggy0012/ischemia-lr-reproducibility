"""Export original QC UMI matrices, invoke official R DecontX, and audit portable output."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,gzip,shutil,gc,re
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS']:
    os.environ[v]='1'
import numpy as np,pandas as pd
from scipy import sparse
from revision_data import load_sample,samples
from repro_decontx_io import write_verified
ROOT=Path(__file__).resolve().parent
SOURCE_BYTES={name:(ROOT/name).read_bytes() for name in ['repro_decontx.py','repro_decontx_run.R','repro_decontx_io.py','repro_decontx_plan.md']}
SOURCE_HASHES={name:hashlib.sha256(value).hexdigest() for name,value in SOURCE_BYTES.items()}
OUT=Path(os.environ.get('DECONTX_OUTPUT_ROOT',ROOT/'repro_decontx'))
RSCRIPT=Path(os.environ.get('DECONTX_RSCRIPT',ROOT/'r_runtime'/'bin'/'Rscript.exe'))
RLIBRARY=Path(os.environ.get('DECONTX_RLIBRARY',ROOT/'repro_decontx'/'rlibrary'))

def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def run(acc,sample):
    dest=OUT/acc/sample;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'audit.json').exists():
        print('CACHED',acc,sample,flush=True);return
    if not (dest/'r_audit.json').exists():
        if not (dest/'input_audit.json').exists():
            (dest/'executed_decontx_run.R').write_bytes(SOURCE_BYTES['repro_decontx_run.R'])
            (dest/'python_wrapper_source_snapshot.py').write_bytes(SOURCE_BYTES['repro_decontx.py'])
            (dest/'io_source_snapshot.py').write_bytes(SOURCE_BYTES['repro_decontx_io.py'])
            x,q,genes=load_sample(acc,sample)
            if not x.has_sorted_indices:x=x.sorted_indices()
            idpath=ROOT/'revision_results'/'identity'/acc/f'{sample}_identity.tsv.gz'
            identity=pd.read_csv(idpath,sep='\t')
            assert np.array_equal(q.barcode,identity.barcode)
            assert q.barcode.is_unique and genes.is_unique
            q=q.copy();q['whole_brain_broad']=identity.whole_brain_broad
            assert q.whole_brain_broad.notna().all()
            q.to_csv(dest/'cells.tsv.gz',sep='\t',index=False)
            pd.DataFrame({'symbol':genes}).to_csv(dest/'genes.tsv.gz',sep='\t',index=False)
            x.data.astype('<i4',copy=False).tofile(dest/'input_data.bin')
            x.indices.astype('<i4',copy=False).tofile(dest/'input_indices.bin')
            x.indptr.astype('<i4',copy=False).tofile(dest/'input_indptr.bin')
            meta={'dataset':acc,'sample':sample,'n_genes':x.shape[1],'n_cells':x.shape[0],
                  'nnz':x.nnz,'raw_total':int(x.sum()),'input_orientation':'genes x cells CSC (same slots as cells x genes CSR)',
                  'identity_sha256':digest(idpath),'label_counts':q.whole_brain_broad.value_counts().to_dict(),
                  'script_sha256':SOURCE_HASHES,
                  'raw_matrix_cache':str(ROOT/'revision_cache'/acc/sample),
                  'input_sha256':{p.name:digest(p) for p in dest.glob('input_*.bin')}}
            (dest/'input_audit.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
            assert digest(dest/'executed_decontx_run.R')==meta['script_sha256']['repro_decontx_run.R']
            assert digest(dest/'python_wrapper_source_snapshot.py')==meta['script_sha256']['repro_decontx.py']
            (dest/'source_snapshot_provenance.json').write_text(json.dumps({
                'snapshot':'executed_decontx_run.R','sha256':digest(dest/'executed_decontx_run.R'),
                'recorded_before_model_execution':True,'matches_sha256_recorded_before_model_execution':True,
                'execution':'Rscript executes this immutable per-sample snapshot, not the mutable source file.'},indent=2),encoding='utf-8')
            (dest/'python_source_snapshot_provenance.json').write_text(json.dumps({
                'snapshot':'python_wrapper_source_snapshot.py','sha256':SOURCE_HASHES['repro_decontx.py'],
                'snapshot_recorded_before_model_execution':True,'matches_sha256_recorded_before_model_execution':True,
                'note':'Source bytes frozen when this Python process started; R executes its per-library snapshot.'},indent=2),encoding='utf-8')
            print('EXPORTED',acc,sample,meta['n_cells'],flush=True)
            del x,q,genes,identity;gc.collect()
        saved_meta=json.loads((dest/'input_audit.json').read_text())
        assert digest(dest/'executed_decontx_run.R')==saved_meta['script_sha256']['repro_decontx_run.R']
        command=[str(RSCRIPT),'--vanilla',str(dest/'executed_decontx_run.R'),str(RLIBRARY),str(dest)]
        with open(dest/'run.log','w',encoding='utf-8') as log:
            proc=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
        if proc.returncode:
            print((dest/'run.log').read_text(encoding='utf-8',errors='replace')[-5000:],flush=True)
            raise RuntimeError(f'R DecontX failed {acc}/{sample}: {proc.returncode}')
    meta=json.loads((dest/'input_audit.json').read_text());audit=json.loads((dest/'r_audit.json').read_text())
    convergence_logs=re.findall(r'Completed iteration:\s*(\d+)\s*\| converge:\s*([0-9.eE+-]+)',(dest/'run.log').read_text(encoding='utf-8',errors='replace'))
    audit['last_logged_iteration']=int(convergence_logs[-1][0]) if convergence_logs else None
    audit['last_logged_max_divergence']=float(convergence_logs[-1][1]) if convergence_logs else None
    audit['convergence_threshold_reached']=bool(convergence_logs and float(convergence_logs[-1][1])<.001)
    vals=np.memmap(dest/'corrected_data.bin',dtype='<f8',mode='r')
    inds=np.memmap(dest/'corrected_indices.bin',dtype='<i4',mode='r')
    ptr=np.memmap(dest/'corrected_indptr.bin',dtype='<i4',mode='r')
    x=sparse.csr_matrix((vals,inds,ptr),shape=(meta['n_cells'],meta['n_genes']),copy=False)
    cont=pd.read_csv(dest/'contamination.tsv.gz',sep='\t');q=pd.read_csv(dest/'cells.tsv.gz',sep='\t')
    assert np.array_equal(cont.barcode,q.barcode)
    assert np.allclose(np.asarray(x.sum(axis=1)).ravel(),cont.corrected_total,rtol=1e-12,atol=1e-7)
    assert np.isfinite(vals).all() and (vals>=0).all()
    (dest/'postprocessing_source_snapshot.py').write_bytes(SOURCE_BYTES['repro_decontx.py'])
    (dest/'io_source_snapshot.py').write_bytes(SOURCE_BYTES['repro_decontx_io.py'])
    io_audit=write_verified(dest,x)
    audit.update({'matrix_npz_orientation':'cells x genes CSR','matrix_market_orientation':'genes x cells',
                  'io_verification':io_audit,
                  'postprocessing_script_sha256':{name:SOURCE_HASHES[name] for name in ['repro_decontx.py','repro_decontx_io.py']},
                  'counts_type':'fractional computationally decontaminated estimates; NOT observed UMI counts',
                  'input':meta,'output_sha256':{name:digest(dest/name) for name in ['corrected_counts.npz','corrected_counts.mtx.gz','genes.tsv.gz','cells.tsv.gz','contamination.tsv.gz','model_estimates.rds']},
                  'python_verification':{'barcode_order':True,'corrected_totals':True,'finite_nonnegative':True}})
    (dest/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    del x,vals,inds,ptr;gc.collect()
    # Remove only task-local redundant transport files, after successful audited outputs.
    for name in ['input_data.bin','input_indices.bin','input_indptr.bin','corrected_data.bin','corrected_indices.bin','corrected_indptr.bin','corrected_counts.mtx']:
        (dest/name).unlink()
    print('VERIFIED',acc,sample,'contamination',audit['contamination_quantiles'],flush=True)

def summarize():
    rows=[];groups=[]
    for path in sorted(OUT.glob('GSE*/GSM*/audit.json')):
        a=json.loads(path.read_text());d=path.parent
        q=pd.read_csv(d/'cells.tsv.gz',sep='\t')
        c=pd.read_csv(d/'contamination.tsv.gz',sep='\t')
        assert np.array_equal(q.barcode,c.barcode)
        c['cell_type']=q.cell_type
        base={'dataset':a['dataset'],'sample':a['sample']}
        def metrics(v):
            return {'n_cells':len(v),'mean_estimated_contamination':v.contamination.mean(),
                    'median_estimated_contamination':v.contamination.median(),
                    'q25_estimated_contamination':v.contamination.quantile(.25),
                    'q75_estimated_contamination':v.contamination.quantile(.75),
                    'q95_estimated_contamination':v.contamination.quantile(.95),
                    'raw_total':v.raw_total.sum(),'corrected_total':v.corrected_total.sum(),
                    'removed_total_fraction':1-v.corrected_total.sum()/v.raw_total.sum()}
        rows.append(base|metrics(c)|{'decontX_version':a['decontX'],'elapsed_seconds':a['elapsed_seconds']})
        for label_column in ['whole_brain_broad','cell_type']:
            for label,v in c.groupby(label_column,observed=True):
                groups.append(base|{'label_source':label_column,'label':label}|metrics(v))
    if rows:
        pd.DataFrame(rows).to_csv(OUT/'sample_contamination_summary.tsv',sep='\t',index=False)
        pd.DataFrame(groups).to_csv(OUT/'cell_type_contamination_summary.tsv',sep='\t',index=False)

if __name__=='__main__':
    for acc in sys.argv[1:] or ['GSE174574','GSE245386']:
        for sample in samples(acc):run(acc,sample)
    summarize()
