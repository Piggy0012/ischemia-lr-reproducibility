"""Recover exact pre-parameterization source snapshots and verify the pre-recorded hashes."""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parent
base=ROOT/'repro_decontx'
source=(ROOT/'repro_decontx_run.R').read_text(encoding='utf-8')
old=''.join(line for line in source.splitlines(keepends=True)
            if not line.startswith('max_iter <- ') and not line.startswith('stopifnot(max_iter %in%'))
old=old.replace('maxIter=max_iter','maxIter=500')
raw=old.encode('utf-8');sha=hashlib.sha256(raw).hexdigest()
out=base/'installers'/'executed_decontx_run_original_500.R';out.write_bytes(raw)
for path in base.glob('GSE*/GSM*/input_audit.json'):
    a=json.loads(path.read_text())
    expected=a['script_sha256']['repro_decontx_run.R']
    assert sha==expected,(path,sha,expected)
    dest=path.parent/'executed_decontx_run.R'
    dest.write_bytes(raw)
    assert hashlib.sha256(dest.read_bytes()).hexdigest()==expected
    (path.parent/'source_snapshot_provenance.json').write_text(json.dumps({
        'snapshot':'executed_decontx_run.R','sha256':sha,
        'recorded_before_model_execution':False,
        'matches_sha256_recorded_before_model_execution':True,
        'reconstruction':'Exact original text recovered by reversing the maxIter parameterization addition; cryptographic hash matches pre-run input audit.'},indent=2),encoding='utf-8')
print('EXACT_PRE_RUN_R_SOURCE_RECOVERED',sha,flush=True)
