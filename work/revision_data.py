"""Memory-mapped per-library raw-count access; all original QC cells retained."""
from pathlib import Path
import sys,json,gc
import numpy as np,pandas as pd
from scipy import sparse
from stream_annotate import chunks

ROOT=Path(__file__).resolve().parent

def cache_sample(acc,sample):
    dest=ROOT/'revision_cache'/acc/sample
    manifest=dest/'complete.json'
    if manifest.exists():return dest
    dest.mkdir(parents=True,exist_ok=True)
    processed=ROOT/'processed'/acc
    q=pd.read_csv(processed/f'{sample}_cells.tsv.gz',sep='\t')
    features=pd.read_csv(processed/f'{sample}_features.tsv.gz',sep='\t')
    mask=q.qc_pass.to_numpy(bool);nc=int(mask.sum());nnz=int(q.loc[mask,'n_genes'].sum())
    xdata=np.lib.format.open_memmap(dest/'data.npy',mode='w+',dtype=np.int32,shape=(nnz,))
    xindices=np.lib.format.open_memmap(dest/'indices.npy',mode='w+',dtype=np.int32,shape=(nnz,))
    pos=0;last=-1
    raw=next((ROOT/'data'/acc/'raw').glob(sample+'*matrix.mtx.gz'))
    it=chunks(raw);ng,oldnc,original_nnz=next(it)
    assert oldnc==len(q) and ng==len(features)
    print('CACHE',acc,sample,'QC',nc,'NNZ',nnz,flush=True)
    for a in it:
        assert (np.diff(a[:,1])>=0).all() and a[0,1]>=last,'Matrix must be cell-column ordered'
        last=a[-1,1];a=a[mask[a[:,1]]];n=len(a)
        xdata[pos:pos+n]=a[:,2];xindices[pos:pos+n]=a[:,0];pos+=n
    assert pos==nnz
    xdata.flush();xindices.flush()
    ptr=np.r_[0,np.cumsum(q.loc[mask,'n_genes'].to_numpy(),dtype=np.int64)]
    np.save(dest/'indptr.npy',ptr.astype(np.int32))
    q=q.loc[mask].copy().reset_index(drop=True)
    q.to_csv(dest/'cells.tsv.gz',sep='\t',index=False)
    features.to_csv(dest/'genes.tsv.gz',sep='\t',index=False)
    meta={'dataset':acc,'sample':sample,'shape':[nc,ng],'nnz':nnz,'source':raw.name,'duplicate_symbols':int(features.symbol.duplicated().sum()),'normalization_total':'all-gene UMI from the original matrix'}
    manifest.write_text(json.dumps(meta,indent=2),encoding='utf-8')
    del xdata,xindices;gc.collect()
    return dest

def load_sample(acc,sample):
    dest=cache_sample(acc,sample);meta=json.loads((dest/'complete.json').read_text())
    x=sparse.csr_matrix((np.load(dest/'data.npy',mmap_mode='r'),np.load(dest/'indices.npy',mmap_mode='r'),np.load(dest/'indptr.npy',mmap_mode='r')),shape=meta['shape'],copy=False)
    q=pd.read_csv(dest/'cells.tsv.gz',sep='\t');genes=pd.read_csv(dest/'genes.tsv.gz',sep='\t').symbol.astype(str)
    if genes.duplicated().any():
        unique=pd.Index(genes.unique());codes=unique.get_indexer(genes)
        collapse=sparse.csr_matrix((np.ones(len(genes),np.int32),(np.arange(len(genes)),codes)),shape=(len(genes),len(unique)))
        x=(x@collapse).tocsr();genes=pd.Series(unique)
    assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(),q.n_umis.to_numpy())
    return x,q,pd.Index(genes)

def samples(acc):
    manifest=pd.read_csv(ROOT/'revision_sample_manifest.tsv',sep='\t')
    selected=manifest.loc[manifest.dataset==acc,'sample'].tolist()
    assert selected and len(selected)==len(set(selected)), 'Missing or duplicate sample manifest'
    return sorted(selected)

if __name__=='__main__':
    for acc in sys.argv[1:] or ['GSE174574','GSE245386']:
        for sample in samples(acc):cache_sample(acc,sample)
