"""One explicit offline smoke test in a new disposable library; no model fit.

Requires the original R 4.6.1 baseline and a previously built environment archive.
The inherited baseline is read-only. All writes stay in the test directory or this
metadata directory. Repeated runs refuse to overwrite an existing test directory.
"""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,zipfile
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
TEST=ROOT/'work/repro_v3_cellchat_environment_restore_test'
archive=ROOT/'work/repro_v3_release_assets/cellchat_environment_windows_R461.zip'
integrity=json.loads(archive.with_name('cellchat_environment_archive_integrity.json').read_text())
assert hashlib.sha256(archive.read_bytes()).hexdigest()==integrity['sha256']
if TEST.exists():raise RuntimeError('Refusing to overwrite an existing restore test directory.')
bundle=TEST/'bundle';bundle.mkdir(parents=True)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for item in z.infolist():assert (bundle/item.filename).resolve().is_relative_to(bundle.resolve())
    z.extractall(bundle)
commands=[
 [sys.executable,str(bundle/'environment/verify_bundle.py'),str(bundle)],
 [str(ROOT/'work/r_runtime/bin/Rscript.exe'),'--vanilla',str(bundle/'environment/restore_offline_overlay.R'),
  str(bundle),str(ROOT/'work/repro_decontx/rlibrary'),str(TEST/'fresh_rlibrary')]]
executed=[]
for i,cmd in enumerate(commands):
    p=subprocess.run(cmd,cwd=ROOT,text=True,encoding='utf-8',errors='replace',capture_output=True)
    (TEST/f'step{i+1}.log').write_text(p.stdout+p.stderr,encoding='utf-8')
    executed.append({'command':cmd,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
    if p.returncode:raise RuntimeError(f'Offline smoke step {i+1} failed; inspect {TEST}/step{i+1}.log')
audit=TEST/'fresh_rlibrary_audit.json';session=TEST/'fresh_rlibrary_sessionInfo.txt'
shutil.copyfile(audit,HERE/'offline_restore_smoke_audit.json')
shutil.copyfile(session,HERE/'offline_restore_smoke_sessionInfo.txt')
(HERE/'offline_restore_execution.json').write_text(json.dumps({
 'status':'passed','archive_tested_sha256':integrity['sha256'],
 'test_directory':str(TEST),'new_library':str(TEST/'fresh_rlibrary'),
 'commands':executed,'scientific_models_fitted':False,
 'archive_metadata_may_be_rebuilt_to_add_this_audit':True},indent=2)+'\n',encoding='utf-8')
print(audit.read_text())
