"""Full native CellChat controlled-resource runner; one library and one R process at a time."""
from pathlib import Path
import argparse, gc, gzip, hashlib, json, os, subprocess, time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[k]='1'
import numpy as np
import pandas as pd
import psutil
from scipy import io, sparse
from revision_data import load_sample, samples

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'repro_v3_cellchat'
PROBE=ROOT/'repro_v3_cellchat_probe'
LIMIT=2*1024**3
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False),encoding='utf-8')
def resource():
    d=pd.read_csv(PROBE/'CellChatDB.mouse_protein_canonical_overlap.tsv',sep='\t',keep_default_na=False)
    d=d[(d.n_rows_per_canonical_pair==1)&d.in_liana_mouseconsensus&d.covered_all_11_LR_cofactor_genes].copy()
    assert len(d)==1548 and not d.duplicated(['canonical_ligand','canonical_receptor']).any()
    p=OUT/'controlled_resource.tsv'
    if p.exists():assert p.read_bytes()==d.to_csv(sep='\t',index=False).encode('utf-8')
    else:d.to_csv(p,sep='\t',index=False)
    return d
def signal_genes(d):
    cx=pd.read_csv(PROBE/'CellChatDB.mouse_complex.tsv',sep='\t',keep_default_na=False).set_index('resource_row_id')
    cf=pd.read_csv(PROBE/'CellChatDB.mouse_cofactor.tsv',sep='\t',keep_default_na=False).set_index('resource_row_id')
    def parts(x,tab):
        if not x:return []
        return [str(v) for v in tab.loc[x].to_numpy() if str(v)] if x in tab.index else [x]
    genes=set()
    for _,r in d.iterrows():
        for c in ['ligand','receptor']:genes.update(parts(r[c],cx))
        for c in ['agonist','antagonist','co_A_receptor','co_I_receptor']:genes.update(parts(r[c],cf))
    return sorted(genes)
def export(acc,sample):
    dest=OUT/'controlled'/acc/sample;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'input_audit.json').exists():return dest
    lr=resource();wanted=signal_genes(lr)
    x,q,genes=load_sample(acc,sample)
    orig=json.loads((ROOT/f'repro_liana_rank_diagnostic/raw/{acc}/primary__{sample}.json').read_text())
    keep=(q.cell_type!='Unassigned').to_numpy()
    counts=q.loc[keep,'cell_type'].value_counts();keep &= q.cell_type.isin(counts[counts>=30].index).to_numpy()
    assert q.loc[keep,'cell_type'].value_counts().to_dict()==orig['context_counts']
    idx=genes.get_indexer(wanted);assert np.all(idx>=0)
    y=x[:,idx][keep].astype(np.float32)
    totals=q.loc[keep,'n_umis'].to_numpy(float)
    y=y.multiply((1e4/totals).astype(np.float32)[:,None]).tocsr();y.data=np.log1p(y.data);y.eliminate_zeros()
    assert np.isfinite(y.data).all() and np.all(y.data>0)
    # Promote exactly represented float32 values for R double storage, avoiding re-normalization.
    y=y.astype(np.float64).T.tocsc()
    meta=q.loc[keep,['barcode','cell_type','n_umis']].copy();meta['sample']=sample
    meta.to_csv(dest/'cells.tsv.gz',sep='\t',index=False)
    pd.DataFrame({'symbol':wanted}).to_csv(dest/'genes.tsv',sep='\t',index=False)
    with gzip.open(dest/'logcp10k_signal.mtx.gz','wb') as handle:io.mmwrite(handle,y,precision=17)
    audit={'dataset':acc,'sample':sample,'condition':orig['condition'],'config':'raw_primary',
      'shape_gene_by_cell':list(y.shape),'nnz':int(y.nnz),'context_counts':orig['context_counts'],
      'n_input_all_genes':len(genes),'normalization':'Original LIANA float32 raw counts -> all-gene total scale 10000 -> log1p, then lossless float64 promotion for R',
      'normalization_total_min':float(totals.min()),'normalization_total_max':float(totals.max()),
      'source_raw_cache_meta_sha256':sha(ROOT/f'revision_cache/{acc}/{sample}/complete.json'),
      'source_primary_audit_sha256':sha(ROOT/f'repro_liana_rank_diagnostic/raw/{acc}/primary__{sample}.json'),
      'resource_sha256':sha(OUT/'controlled_resource.tsv'),'frozen_plan_sha256':sha(ROOT/'repro_v3_independent_framework_plan.md'),
      'files_sha256':{name:sha(dest/name) for name in ['cells.tsv.gz','genes.tsv','logcp10k_signal.mtx.gz']}}
    dump(dest/'input_audit.json',audit)
    del x,q,y;gc.collect()
    print('EXPORTED',sample,audit['shape_gene_by_cell'],audit['nnz'],flush=True)
    return dest
def run(acc,sample):
    assert json.loads((OUT/'installed_package_verified.json').read_text())['status']=='complete_native_package_verified'
    dest=export(acc,sample)
    if (dest/'run_telemetry.json').exists():
        prior=json.loads((dest/'run_telemetry.json').read_text())
        if prior['status']=='complete':print('ALREADY COMPLETE',sample,flush=True);return prior
    cmd=[str(ROOT/'r_runtime/bin/Rscript.exe'),'--vanilla',str(ROOT/'repro_v3_cellchat_run.R'),str(dest)]
    t0=time.monotonic();peak=0;trace=[];reason=None
    with (dest/'model.log').open('w',encoding='utf-8') as log:
        p=subprocess.Popen(cmd,cwd=ROOT.parent,stdout=log,stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                proc=psutil.Process(p.pid);children=proc.children(recursive=True)
                rss=sum(c.memory_info().rss for c in [proc]+children if c.is_running());peak=max(peak,rss)
                trace.append({'elapsed_seconds':round(time.monotonic()-t0,2),'rss_bytes':rss})
                if rss>LIMIT:
                    reason='Exceeded 2 GiB actual process-tree RSS';
                    for c in children:c.kill()
                    p.kill();break
            except psutil.NoSuchProcess:pass
            time.sleep(1)
        rc=p.wait()
    rep={'status':'complete' if rc==0 and reason is None else 'failed','dataset':acc,'sample':sample,
         'returncode':rc,'stop_reason':reason,'peak_process_tree_rss_bytes':peak,'elapsed_wall_seconds':time.monotonic()-t0,
         'limit_bytes':LIMIT,'R_script_sha256':sha(ROOT/'repro_v3_cellchat_run.R'),'runner_sha256':sha(__file__),
         'plan_sha256':sha(ROOT/'repro_v3_independent_framework_plan.md'),'installed_package_audit_sha256':sha(OUT/'installed_package_verified.json')}
    pd.DataFrame(trace).to_csv(dest/'memory_trace.tsv',sep='\t',index=False)
    if rep['status']=='complete':rep['output_sha256']={n:sha(dest/n) for n in ['model_audit.json','full_network.tsv.gz','native_network.rds']}
    dump(dest/'run_telemetry.json',rep);print(json.dumps(rep),flush=True)
    if rep['status']!='complete':raise RuntimeError(f'Native CellChat failed: {dest}/model.log')
    return rep
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--dataset');ap.add_argument('--sample');ap.add_argument('--all',action='store_true');ap.add_argument('--export-only',action='store_true');args=ap.parse_args()
    if args.all:
        for acc in ['GSE174574','GSE245386']:
            for s in samples(acc):run(acc,s)
    else:
        assert args.dataset and args.sample
        (export if args.export_only else run)(args.dataset,args.sample)
