"""Bind exact historical Python sources to already recorded pre-model SHA256 values."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'repro_decontx'
current=(ROOT/'repro_decontx.py').read_text(encoding='utf-8')
middle=current.replace('from repro_decontx_io import write_verified\n','')
middle=middle.replace('    io_audit=write_verified(dest,x)\n',
    "    sparse.save_npz(dest/'corrected_counts.npz',x,compressed=True)\n"
    "    with open(dest/'corrected_counts.mtx','rb') as src,gzip.open(dest/'corrected_counts.mtx.gz','wb',compresslevel=5) as dst:\n"
    "        shutil.copyfileobj(src,dst,8*1024*1024)\n")
middle=''.join(line for line in middle.splitlines(keepends=True)
    if "'io_verification':io_audit" not in line and "'postprocessing_script_sha256'" not in line)
old=(BASE/'installers'/'executed_python_original.py').read_bytes()
candidates={hashlib.sha256(raw).hexdigest():raw for raw in [old,middle.encode('utf-8')]}
folder=BASE/'source_snapshots';folder.mkdir(exist_ok=True)
records=[]
for path in BASE.glob('GSE*/GSM*/input_audit.json'):
    a=json.loads(path.read_text());sha=a['script_sha256']['repro_decontx.py']
    assert sha in candidates,(path,sha,list(candidates))
    raw=candidates[sha]
    assert hashlib.sha256(raw).hexdigest()==sha
    common=folder/f'repro_decontx_{sha}.py';common.write_bytes(raw)
    dest=path.parent/'python_wrapper_source_snapshot.py';dest.write_bytes(raw)
    provenance={'snapshot':'python_wrapper_source_snapshot.py','sha256':sha,
        'matches_sha256_recorded_before_model_execution':True,
        'snapshot_recorded_before_model_execution':False,
        'source_origin':'Preserved original Python source' if raw==old else 'Exact immutable-R-snapshot version recovered by reversing the later I/O-only patch',
        'note':'Code content snapshot; the original invocation used the source file in the work directory. Postprocessing code changes are tracked separately.'}
    (path.parent/'python_source_snapshot_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    records.append({'dataset':a['dataset'],'sample':a['sample'],'sha256':sha,'snapshot':str(common.relative_to(BASE))})
(folder/'historical_python_source_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print('VERIFIED_HISTORICAL_PYTHON_SNAPSHOTS',len(records),list(candidates),flush=True)
