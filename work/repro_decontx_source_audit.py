"""Verify model-launch and output-transport source hashes without loading count matrices."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

records=[]
for directory in [ROOT/'repro_decontx',ROOT/'repro_decontx_maxiter2000']:
    if not directory.exists():continue
    for p in directory.glob('GSE*/GSM*/audit.json'):
        a=json.loads(p.read_text());dest=p.parent
        mapping={'repro_decontx.py':dest/'python_wrapper_source_snapshot.py',
                 'repro_decontx_run.R':dest/'executed_decontx_run.R',
                 'repro_decontx_io.py':dest/'io_source_snapshot.py',
                 'repro_decontx_plan.md':ROOT/'repro_decontx_plan.md'}
        for name,expected in a['input']['script_sha256'].items():
            assert digest(mapping[name])==expected,(p,name,'model-launch source differs')
        post_mapping={'repro_decontx.py':dest/'postprocessing_source_snapshot.py',
                      'repro_decontx_io.py':dest/'io_source_snapshot.py'}
        for name,expected in a.get('postprocessing_script_sha256',{}).items():
            assert digest(post_mapping[name])==expected,(p,name,'postprocessing source differs')
        records.append({'directory':directory.name,'dataset':a['dataset'],'sample':a['sample'],
                        'maxIter':a['requested_parameters']['maxIter'],
                        'model_sources_verified':len(a['input']['script_sha256']),
                        'postprocessing_sources_verified':len(a.get('postprocessing_script_sha256',{})),
                        'audit_recovery_present':'audit_recovery' in a})
report={'all_declared_source_hashes_match':True,'completed_outputs_checked':len(records),'records':records}
(ROOT/'repro_decontx'/'source_integrity_check.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
