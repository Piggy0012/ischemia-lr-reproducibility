"""Collect compact relocated-run evidence; preserve local before-output backups."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import shutil
import pandas as pd

WORK = Path(__file__).resolve().parent
BASE = WORK.parent
STAGE = WORK / 'repro_delivery/analysis'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def differences(a, b, path=''):
    if isinstance(a, dict) and isinstance(b, dict):
        answer = []
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                answer.append({'field': path + '/' + k, 'before': a.get(k), 'after': b.get(k)})
            else:
                answer.extend(differences(a[k], b[k], path + '/' + k))
        return answer
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in differences(x, y, path + '/' + str(i))]
    return [] if a == b else [{'field': path, 'before': a, 'after': b}]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_id')
    ap.add_argument('--failed-path-diagnostic', action='store_true')
    ap.add_argument('--after-root', type=Path, help='An immutable next-run before/ snapshot for an earlier diagnostic.')
    args = ap.parse_args()
    runs = STAGE / 'work/repro_cached_runs'
    run = (runs / args.run_id).resolve()
    assert run.parent == runs.resolve() and run.is_dir()
    audit = json.loads((run / 'cached_pipeline_audit.json').read_text())
    after_root = args.after_root.resolve() if args.after_root else STAGE
    assert after_root.is_relative_to(STAGE.resolve())
    if not args.failed_path_diagnostic:
        assert audit['status'] == 'passed' and audit['exit_code'] == 0
        assert len(audit['stages']) == 12 and not audit['pending_stages']
        assert not audit['scientific_table_changes'] and not audit['protected_input_changes']
    name = 'relocated_path_failure' if args.failed_path_diagnostic else 'relocated_cached_pipeline'
    destination = STAGE / 'outputs/reproducibility_v2/validation' / name
    destination.mkdir(parents=True, exist_ok=True)
    for p in run.iterdir():
        if p.is_file() and p.suffix in {'.json', '.tsv', '.log'}:
            shutil.copy2(p, destination / p.name)
    shutil.copytree(run / 'validation_after', destination / 'validation_after', dirs_exist_ok=True)
    comparisons = pd.read_csv(run / 'semantic_comparison.tsv', sep='\t')
    details = []
    for r in comparisons[comparisons.status.eq('audit_metadata_changed')].itertuples(index=False):
        # The source of the after metadata is the relocated stage, never root.
        rel = Path(r.path)
        before, after = run / 'before' / rel, after_root / rel
        assert sha(before) == r.before_sha256 and sha(after) == r.after_sha256
        details.append({'path': r.path, 'differences': differences(json.loads(before.read_text()), json.loads(after.read_text()))})
    (destination / 'metadata_differences.json').write_text(json.dumps(details, indent=2, ensure_ascii=False), encoding='utf-8')
    review = ('This initial relocated cached run completed every computation with no scientific table or protected-source change, '
              'but failed post-validation because the validator selected an outer /work/ path marker before the inner /outputs/ marker. '
              'The validator was corrected to choose the innermost archive-root marker; six root, relocated and relative path cases passed. '
              'This directory deliberately preserves the failed audit and is not evidence of a successful release.'
              if args.failed_path_diagnostic else
              'All 12 cached stages passed after relocation, followed by complete validation with no pending stages. '
              'All scientific table schemas, row/column order, missingness and values were preserved; protected input and historical-test hashes were unchanged. '
              'The analysis scripts and all data sources resolved inside the relocated stage. The existing original frozen Python environment was shared; '
              'this is a directory-relocation test, not a fresh-machine environment installation test.')
    (destination / 'README.md').write_text('# Relocated cached reconstruction evidence\n\n' + review + '\n\n'
        'No large cell matrices, model fits, software installation, network access, or figure regeneration occurred. '
        'The complete local before-output duplicate remains in the recorded run directory; it is excluded from the release ZIP. '
        'Compact plans, logs, validations, source hashes, semantic comparison and explicit JSON metadata differences are retained here. '
        'The successful comparison includes floating-point values within the recorded strict tolerance; gzip container timestamps are not scientific changes.\n', encoding='utf-8')
    payload = [{'path': p.relative_to(destination).as_posix(), 'bytes': p.stat().st_size, 'sha256': sha(p)}
               for p in sorted(destination.rglob('*')) if p.is_file() and p.name != 'evidence_manifest.json']
    manifest = {'status': audit['status'], 'collected_utc': datetime.now(timezone.utc).isoformat(),
                'source_run_relative_to_relocated_release_root': run.relative_to(STAGE).as_posix(),
                'source_run_relative_to_original_workspace': run.relative_to(BASE).as_posix(),
                'directory_relocation_test': True, 'fresh_runtime_installation_test': False,
                'after_metadata_source_relative_to_stage': after_root.relative_to(STAGE).as_posix(),
                'semantic_status_counts': comparisons.status.value_counts().to_dict(),
                'copied_artifacts': payload, 'manifest_excludes_itself': True,
                'local_before_backup_included': False, 'collector_script_sha256': sha(Path(__file__))}
    (destination / 'evidence_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    root_destination = BASE / 'outputs/reproducibility_v2/validation' / name
    shutil.copytree(destination, root_destination, dirs_exist_ok=True)
    print(json.dumps({'status': audit['status'], 'evidence': str(root_destination),
                      'file_count': len(payload) + 1, 'semantic_status_counts': manifest['semantic_status_counts']}, indent=2))


if __name__ == '__main__':
    main()
