"""Audit DecontX source coverage in the archive staging rules, without copying data."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
STAGE = ROOT / 'repro_delivery/analysis/work'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    staging = ROOT / 'repro_stage_archive.py'
    text = staging.read_text(encoding='utf-8')
    assert "['.json','.tsv','.md','.R','.py','.txt']" in text
    required = set(ROOT.glob('repro_decontx*.*'))
    for directory in ['repro_decontx', 'repro_decontx_maxiter2000']:
        base = ROOT / directory
        required.update(p for p in base.glob('*') if p.is_file() and p.suffix in ['.json', '.tsv', '.md', '.R', '.py', '.txt'])
        for sample in base.glob('GSE*/GSM*'):
            required.update(p for p in sample.glob('*') if p.is_file() and p.suffix in ['.json', '.R', '.py', '.log', '.txt'])
    required.update(p for p in (ROOT / 'repro_decontx/source_snapshots').rglob('*') if p.is_file())
    required.update(p for p in (ROOT / 'repro_decontx/installers').glob('*')
                    if p.is_file() and (p.suffix in ['.json', '.R', '.py', '.txt'] or p.name.endswith(('_LICENSE', '_DESCRIPTION'))))
    # The final report is generated below and will itself match repro_*.json.
    required.discard(ROOT / 'repro_decontx_archive_source_check.json')
    provenance_path = ROOT / 'repro_decontx/GSE245386/GSM7841720/audit_recovery_provenance.json'
    recovery = json.loads(provenance_path.read_text())
    for name, digest in recovery['recovery_script_sha256'].items():
        assert sha(ROOT / name) == digest
    post = recovery['postprocessing_after_recovery']
    assert sha(provenance_path.parent / post['snapshot']) == post['script_sha256']
    selection = json.loads((ROOT / 'repro_decontx_maxiter2000/selection.json').read_text())
    assert sha(ROOT / 'repro_decontx_maxiter2000/extension_driver_source_snapshot.py') == selection['extension_driver_sha256']
    assert sha(ROOT / 'repro_decontx_convergence_amendment.md') == selection['amendment_sha256']
    assert sha(ROOT / 'repro_decontx_maxiter2000/analysis_plan_snapshot.md') == selection['amendment_sha256']
    records = []
    for p in sorted(required):
        rel = p.relative_to(ROOT)
        staged = STAGE / rel
        records.append({'path': str(rel).replace('\\', '/'), 'sha256': sha(p),
                        'staged_exists': staged.exists(),
                        'staged_hash_matches': sha(staged) == sha(p) if staged.exists() else None})
    report = {
        'status': 'copy_rules_checked', 'staging_script_sha256': sha(staging),
        'stage_exists': STAGE.exists(),
        'copy_rules_cover_required_decontx_sources': True,
        'extension_root_py_snapshot_rule_present': True,
        'recovery_and_extension_additional_hashes_match': True,
        'required_source_and_audit_files': len(records),
        'all_existing_staged_files_match': all(r['staged_hash_matches'] for r in records if r['staged_exists']) if any(r['staged_exists'] for r in records) else None,
        'actual_archive_validation_still_required': True,
        'note': 'Static rules and local sources were inspected; this is not a claim that a release archive exists or has been verified.',
        'records': records,
    }
    (ROOT / 'repro_decontx_archive_source_check.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'records'}, indent=2))


if __name__ == '__main__':
    main()
