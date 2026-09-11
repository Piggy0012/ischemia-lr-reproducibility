"""Index the captured effective R libraries and optionally archive the exact 42 ZIPs.

No package installation, network request or model fitting. Metadata is suitable
for Git; binary assets are written outside the repository as a separate archive.
"""
from pathlib import Path
import argparse, csv, hashlib, json, shutil, zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CELL = ROOT / 'work/repro_v3_cellchat'
BASE = ROOT / 'work/repro_v3_environment'
ASSETS = ROOT / 'work/repro_v3_release_assets'

def digest(path, algo='sha256'):
    h = hashlib.new(algo)
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def save(path, value): Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
def table(path):
    with Path(path).open(encoding='utf-8', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
def tsv(path, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)

def main(archive=False):
    cap=load(HERE/'capture_audit.json'); desc=load(HERE/'R_package_descriptions.json')
    lock=load(HERE/'renv.lock'); eff={r['Package']:r for r in table(HERE/'R_packages_analysis_effective.tsv')}
    prior=load(BASE/'R_package_descriptions.json')
    source=load(CELL/'download_manifest.json')
    assert len(eff)==287 and len(lock['Packages'])==274 and len(desc)==288
    assert cap['overlay_existing_packages_with_changed_versions']==[]
    provenance=HERE/'provenance';provenance.mkdir(exist_ok=True)
    for name in ['download_manifest.json','installed_package_verified.json','environment_summary.json','analysis_sessionInfo.txt']:
        shutil.copyfile(CELL/name,provenance/name)
    save(provenance/'source_binding.json',{'source_files_sha256':{
        p.relative_to(ROOT).as_posix():digest(p) for p in [
          ROOT/'work/repro_v3_cellchat_run.py',ROOT/'work/repro_v3_cellchat_run.R',
          ROOT/'work/repro_v3_cellchat_install_verify.R',ROOT/'work/repro_v3_cellchat_capture_environment.R',
          ROOT/'work/repro_v3_independent_framework_plan.md',BASE/'renv.lock']},
        'native_CellChat_commit':'75253cd0c9e68410e6e721a6d3a0419a1d7e358f',
        'native_binary_provider':'https://blaserlab.r-universe.dev',
        'official_source_repository':'https://github.com/jinworks/CellChat',
        'package_source_verification':'See installed_package_verified.json: 47 function bodies/formals and entire mouse DB; not whole-package byte identity to a source build.'})
    bins=[]
    for r in source['packages']:
        p=ROOT/Path(r['local_path']); d=desc[r['package']]
        assert p.is_file() and p.stat().st_size==r['bytes'] and digest(p)==r['sha256']
        assert d['Version']==r['version'] and zipfile.is_zipfile(p)
        with zipfile.ZipFile(p) as z: assert z.testzip() is None
        bins.append({'Package':r['package'],'Version':r['version'],
          'archive_path':'binaries/'+p.name,'original_url':r['url'],'bytes':r['bytes'],
          'sha256':r['sha256'],'md5':digest(p,'md5')})
    assert len(bins)==42 and {r['Package'] for r in bins}==set(cap['overlay_new_packages'])
    tsv(HERE/'offline_binary_manifest.tsv',bins)
    rows=[]; unchanged=[]
    for p,d in sorted(desc.items()):
        scope='capture_tool_only' if p=='renv' else ('cellchat_overlay' if p in cap['overlay_new_packages'] else 'baseline_or_R_runtime')
        if scope=='baseline_or_R_runtime':
            assert prior[p]['Version']==d['Version']
            assert digest(HERE/'R_DESCRIPTION'/f'{p}.DESCRIPTION')==digest(BASE/'R_DESCRIPTION'/f'{p}.DESCRIPTION')
            unchanged.append(p)
        if p in lock['Packages']: assert lock['Packages'][p]['Version']==d['Version']
        rows.append({'Package':p,'Version':d['Version'],'scope':scope,'Priority':d.get('Priority',''),
          'Built':d.get('Built',''),'Repository':d.get('Repository',''),'RemoteUrl':d.get('RemoteUrl',''),
          'RemoteRef':d.get('RemoteRef',''),'RemoteSha':d.get('RemoteSha',''),
          'description_sha256':digest(HERE/'R_DESCRIPTION'/f'{p}.DESCRIPTION'),
          'in_analysis_libpaths':int(p in eff),'in_renv_lock':int(p in lock['Packages'])})
    tsv(HERE/'package_lock.tsv',rows)
    audit={'status':'captured_lock_and_original_binaries_verified','effective_analysis_packages':287,
      'additional_capture_only_package':'renv 1.2.4','locked_nonbase_packages':274,
      'inherited_packages_description_byte_identical_to_baseline':len(unchanged),
      'new_overlay_packages':42,'baseline_package_versions_changed':0,
      'binary_total_bytes':sum(r['bytes'] for r in bins),'binary_sha256_and_inner_zip_CRC_verified':42,
      'full_environment_restored_from_empty_library':False,
      'metadata_git_suitable':True,'binaries_git_suitable':False,
      'lock_restore_boundary':'renv.lock records installed sources and versions; a full source restore was not tested. Exact archived Windows binaries pin the 42-package overlay, with the inherited baseline separately required.',
      'capture_validation_fix':'Exact DESCRIPTION Version comparison avoids packageVersion hyphen normalization. Quoted TSV fields preserve multiline DESCRIPTION metadata. These capture-only validator/serialization fixes changed no installed package or scientific definition.'}
    smoke=HERE/'offline_restore_smoke_audit.json'
    if smoke.exists():
        s=load(smoke);assert s['status']=='offline_overlay_installed_and_namespace_loaded'
        audit['offline_overlay_restore_tested']=True
        audit['offline_restore_expected_versions_checked']=s['effective_analysis_package_versions_checked']
        audit['offline_restore_scope']='42 archived binary packages freshly installed in an isolated library, inheriting the existing baseline; native namespace loaded, no new model fit.'
    else:audit['offline_overlay_restore_tested']=False
    python=HERE/'python_current_audit.json'
    if python.exists():
        p=load(python);assert p['status']=='current_environment_and_optional_audit_dependencies_verified'
        audit['python_audit_closure_verified']=True
        audit['python_current_pins']=p['current_pins']
        audit['python_pin_delta_from_baseline']=p['added']
        audit['python_independent_review_dependency_closure_packages']=p['independent_review_required_dependency_closure_packages']
    save(HERE/'environment_lock_audit.json',audit)
    excluded={'environment_manifest.json'}
    payload=[p for p in HERE.rglob('*') if p.is_file() and p.name not in excluded and '__pycache__' not in p.parts]
    entries=[{'path':p.relative_to(HERE).as_posix(),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(payload)]
    save(HERE/'environment_manifest.json',{'self_excluded':True,'files':entries,'count':len(entries),
        'bytes':sum(r['bytes'] for r in entries)})
    if archive:
        ASSETS.mkdir(exist_ok=True)
        dest=ASSETS/'cellchat_environment_windows_R461.zip'
        files=[(p,'environment/'+p.relative_to(HERE).as_posix()) for p in sorted(payload+[HERE/'environment_manifest.json'])]
        files += [(ROOT/Path(r['local_path']),'binaries/'+Path(r['local_path']).name) for r in source['packages']]
        with zipfile.ZipFile(dest,'w',allowZip64=True) as z:
            for p,name in files: z.write(p,name,compress_type=zipfile.ZIP_STORED if p.suffix=='.zip' else zipfile.ZIP_DEFLATED,compresslevel=None if p.suffix=='.zip' else 1)
        with zipfile.ZipFile(dest) as z:
            assert z.testzip() is None
            for p,name in files: assert hashlib.sha256(z.read(name)).hexdigest()==digest(p)
        rep={'status':'all_member_sha256_and_CRC_passed','archive':dest.name,'bytes':dest.stat().st_size,
          'sha256':digest(dest),'members':len(files),'binary_packages':42,
          'contains_R_runtime':False,'contains_baseline_binary_packages':False,
          'scope':'Windows R 4.6.1 offline CellChat overlay plus complete effective-library metadata; baseline dependencies are required separately.'}
        save(ASSETS/'cellchat_environment_archive_integrity.json',rep);print(json.dumps(rep,indent=2))
    else: print(json.dumps(audit,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--archive',action='store_true');a=p.parse_args();main(a.archive)
