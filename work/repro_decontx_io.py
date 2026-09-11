"""Lossless portable sparse-matrix transport using low compression and exact readback."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name]='1'
import gzip,shutil,zipfile
from io import BytesIO
from pathlib import Path
import numpy as np
from scipy import sparse,io
import scipy
from threadpoolctl import threadpool_limits,threadpool_info

def save_npz_level1(path,x):
    assert sparse.isspmatrix_csr(x)
    arrays={'indices':x.indices,'indptr':x.indptr,'format':np.asarray(b'csr'),
            'shape':np.asarray(x.shape),'data':x.data}
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as z:
        for name,array in arrays.items():
            with z.open(name+'.npy','w',force_zip64=True) as f:
                np.lib.format.write_array(f,np.asarray(array),allow_pickle=False)

def equal_csr(x,y):
    return (x.shape==y.shape and np.array_equal(x.indptr,y.indptr)
            and np.array_equal(x.indices,y.indices) and np.array_equal(x.data,y.data))

def write_verified(dest,x):
    dest=Path(dest)
    save_npz_level1(dest/'corrected_counts.npz',x)
    check=sparse.load_npz(dest/'corrected_counts.npz')
    assert equal_csr(x,check),'NPZ full matrix readback differs'
    del check
    mtx=dest/'corrected_counts.mtx'
    generated=not mtx.exists()
    if generated:
        # Load the lazily imported I/O backend before applying threadpoolctl.
        io.mminfo(BytesIO(b'%%MatrixMarket matrix coordinate real general\n1 1 0\n'))
        with threadpool_limits(limits=1):
            pools=[p for p in threadpool_info() if p['internal_api']=='scipy_mmio']
            assert pools and all(p['num_threads']==1 for p in pools)
            io.mmwrite(mtx,x.T,field='real',precision=17,symmetry='general')
            check=io.mmread(mtx,spmatrix=True).T.tocsr()
        assert equal_csr(x,check),'17-digit Matrix Market full matrix readback differs'
        del check
    with open(mtx,'rb') as src,gzip.open(dest/'corrected_counts.mtx.gz','wb',compresslevel=1) as dst:
        shutil.copyfileobj(src,dst,8*1024*1024)
    with gzip.open(dest/'corrected_counts.mtx.gz','rt') as f:
        assert next(f).strip()=='%%MatrixMarket matrix coordinate real general'
        line=next(f)
        while line.startswith('%'):line=next(f)
        assert tuple(map(int,line.split()))==(x.shape[1],x.shape[0],x.nnz)
    return {'npz_compression_level':1,'gzip_compression_level':1,'npz_full_readback_exact':True,
            'matrix_market_writer':'scipy.io.mmwrite, precision=17' if generated else 'preserved original R Matrix::writeMM output',
            'matrix_market_full_readback_exact':True if generated else None,
            'matrix_market_gzip_header_and_dimensions_verified':True,
            'io_threads':1,'scipy_version':scipy.__version__}

if __name__=='__main__':
    import tempfile,json
    with tempfile.TemporaryDirectory(prefix='decontx_io_') as tmp:
        data=np.array([np.pi,np.nextafter(1.,2.),1.e-24,999999.9999999999,0.00012345678901234567])
        x=sparse.csr_matrix((data,([0,0,1,2,3],[0,2,1,3,0])),shape=(5,4))
        result=write_verified(tmp,x)
        print(json.dumps(result,indent=2))
