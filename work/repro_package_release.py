"""Package an already-staged, completely validated release. Default: dry run.

Never invokes staging, downloads, model fitting or account/publication tools.
Only the explicit analysis stage and convergence-selected model directories
are eligible inputs. release_integrity.json is external to avoid hash cycles.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

WORK = Path(__file__).resolve().parent
BASE = WORK.parent
DEFAULT_STAGE = WORK / 'repro_delivery/analysis'
DEFAULT_OUTPUT = WORK / 'repro_delivery/release'
DENIED_SEGMENTS = {'.git', '.ssh', '.aws', '.azure', '.config', '.cache', '__pycache__',
                   'venv', '.venv', 'revision_env', 'node_modules', 'rlibrary', 'site-packages',
                   'runtime', 'r-runtime', 'credentials', 'credential'}
DENIED_NAMES = {'.env', '.netrc', '.npmrc', '.pypirc', 'credentials.json', 'credentials.toml',
                'token.json', 'tokens.json', 'auth.json', 'id_rsa', 'id_ed25519', 'cookies.txt'}
STORED_ENDINGS = ('.gz', '.npz', '.png', '.jpg', '.jpeg', '.pdf', '.docx', '.xlsx', '.pptx', '.rds', '.rda', '.zip')
MANIFEST_NAME = 'MANIFEST_SHA256.tsv'
INTEGRITY_NAME = 'release_integrity.json'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inside(path, root):
    return Path(path).resolve().is_relative_to(Path(root).resolve())


def assert_safe(path, root):
    p, root = Path(path), Path(root)
    if p.is_symlink() or not inside(p, root):
        raise ValueError('Input escapes its declared source directory: ' + str(p))
    rel = p.relative_to(root)
    lower = {x.lower() for x in rel.parts}
    if lower & DENIED_SEGMENTS or p.name.lower() in DENIED_NAMES or p.name.lower().startswith('.env.'):
        raise ValueError('Runtime/account file is not an eligible release input: ' + str(rel))
    if p.suffix.lower() in {'.pyc', '.pyo', '.pem', '.key', '.exe', '.dll', '.msi', '.pyd'}:
        raise ValueError('Executable runtime or account-key file is not an eligible release input: ' + str(rel))


def stage_files(stage):
    if not stage.is_dir():
        return []
    result = []
    for p in sorted(stage.rglob('*')):
        if not p.is_file():
            continue
        if p.relative_to(stage).parts[:2] == ('work', 'repro_cached_runs'):
            continue  # Local before-output duplicates; compact audit evidence is in outputs/validation.
        if '__pycache__' in p.parts or p.suffix.lower() in {'.pyc', '.pyo'}:
            continue  # Interpreter caches are reproducible runtime byproducts.
        # Generated manifests are external/inserted as ZIP members, never source
        # payload files. Other archives indicate stale/re-entrant staging.
        if p.name in {MANIFEST_NAME, INTEGRITY_NAME}:
            continue
        if p.suffix.lower() == '.zip':
            raise ValueError('An archive is present inside the analysis stage: ' + str(p))
        assert_safe(p, stage)
        result.append({'path': p, 'member': p.relative_to(stage).as_posix(), 'role': 'code_results'})
    return result


def model_files(strict=False):
    """Read only the 11 named source audits; never discover arbitrary directories."""
    import pandas as pd
    from repro_decontx_source import resolve
    manifest = pd.read_csv(WORK / 'revision_sample_manifest.tsv', sep='\t')
    manifest = manifest[['dataset', 'sample']].drop_duplicates()
    if len(manifest) != 11 or not manifest['sample'].is_unique:
        raise AssertionError('Expected 11 original animal libraries')
    files, selections, pending = {}, [], []
    for row in manifest.itertuples(index=False):
        try:
            selected, decision = resolve(row.dataset, row.sample)
            choices = [('selected', selected)]
            if decision['numerical_extension_selected']:
                choices.append(('original_default_500_nonconvergent', WORK / 'repro_decontx' / row.dataset / row.sample))
            for role, folder in choices:
                allowed = WORK / ('repro_decontx_maxiter2000' if folder.parts[-3] == 'repro_decontx_maxiter2000' else 'repro_decontx')
                if folder.parent.parent not in {WORK / 'repro_decontx', WORK / 'repro_decontx_maxiter2000'}:
                    raise ValueError('Unexpected model source: ' + str(folder))
                audit = json.loads((folder / 'audit.json').read_text(encoding='utf-8'))
                if audit['dataset'] != row.dataset or audit['sample'] != row.sample:
                    raise AssertionError('Model audit identity mismatch')
                if not all(audit['python_verification'].values()):
                    raise AssertionError('Model matrix verification failed')
                for name in ['corrected_counts.npz', 'genes.tsv.gz', 'cells.tsv.gz', 'audit.json']:
                    path = folder / name
                    if not path.is_file():
                        raise FileNotFoundError('Missing required matrix package input: ' + str(path))
                    assert_safe(path, allowed)
                    member = 'work/' + path.relative_to(WORK).as_posix()
                    expected = audit.get('output_sha256', {}).get(name)
                    if name != 'audit.json' and not expected:
                        raise AssertionError('Missing model output hash: ' + member)
                    files[member] = {'path': path, 'member': member, 'role': role,
                                     'expected_sha256': expected,
                                     'dataset': row.dataset, 'sample': row.sample}
                # These trace source order/model settings without duplicating MTX.
                for name in ['input_audit.json', 'r_audit.json', 'sessionInfo.txt', 'installed_packages.tsv',
                             'source_snapshot_provenance.json', 'executed_decontx_run.R',
                             'python_wrapper_source_snapshot.py', 'python_source_snapshot_provenance.json']:
                    path = folder / name
                    if path.is_file():
                        assert_safe(path, allowed)
                        files['work/' + path.relative_to(WORK).as_posix()] = {
                            'path': path, 'member': 'work/' + path.relative_to(WORK).as_posix(), 'role': role,
                            'dataset': row.dataset, 'sample': row.sample}
            selections.append({'dataset': row.dataset, 'sample': row.sample, **decision})
        except (AssertionError, FileNotFoundError) as exc:
            pending.append({'dataset': row.dataset, 'sample': row.sample, 'reason': str(exc)})
    if strict and (pending or len(selections) != 11):
        raise RuntimeError('Corrected matrix package is incomplete: ' + json.dumps(pending))
    return sorted(files.values(), key=lambda x: x['member']), selections, pending


def record_files(files):
    rows = []
    for item in files:
        p = item['path']
        stat = p.stat()
        digest = sha(p)
        if item.get('expected_sha256') and item['expected_sha256'] != digest:
            raise AssertionError('File differs from the completed model audit: ' + item['member'])
        rows.append({'path': item['member'], 'bytes': stat.st_size, 'sha256': digest, 'role': item['role'],
                     'compression': 'stored' if p.name.lower().endswith(STORED_ENDINGS) else 'deflate_level_1'})
    if len({x['path'] for x in rows}) != len(rows):
        raise AssertionError('Duplicate archive member paths')
    return rows


def manifest_bytes(rows):
    text = io.StringIO(newline='')
    writer = csv.DictWriter(text, fieldnames=['path', 'bytes', 'sha256', 'role', 'compression'], delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return text.getvalue().encode('utf-8')


def verify_zip(path, rows, manifest):
    expected = {x['path']: x for x in rows}
    with zipfile.ZipFile(path, 'r') as archive:
        if set(archive.namelist()) != set(expected) | {MANIFEST_NAME}:
            raise AssertionError('ZIP member inventory differs from manifest')
        if len(archive.namelist()) != len(expected) + 1:
            raise AssertionError('Duplicate ZIP member')
        bad = archive.testzip()  # Explicit CRC verification for every member.
        if bad:
            raise AssertionError('ZIP CRC failed: ' + bad)
        if archive.read(MANIFEST_NAME) != manifest:
            raise AssertionError('Embedded manifest differs from external manifest')
        for name, record in expected.items():
            info = archive.getinfo(name)
            if info.file_size != record['bytes']:
                raise AssertionError('ZIP member size differs: ' + name)
            digest = hashlib.sha256()
            with archive.open(name) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    digest.update(block)
            if digest.hexdigest() != record['sha256']:
                raise AssertionError('ZIP member SHA256 differs: ' + name)
    return {'zip_crc_all_members_passed': True, 'zip_member_sha256_all_passed': True,
            'embedded_manifest_matches_external': True, 'manifest_payload_files': len(rows),
            'zip_member_count_including_manifest': len(rows) + 1}


def write_zip(destination, files, rows):
    if destination.exists() or Path(str(destination) + '.partial').exists():
        raise FileExistsError('Release files are immutable; choose a new output directory: ' + str(destination))
    partial = Path(str(destination) + '.partial')
    manifest = manifest_bytes(rows)
    with zipfile.ZipFile(partial, 'x', allowZip64=True) as archive:
        for item, record in zip(files, rows):
            p = item['path']
            compression = zipfile.ZIP_STORED if record['compression'] == 'stored' else zipfile.ZIP_DEFLATED
            archive.write(p, arcname=item['member'], compress_type=compression,
                          compresslevel=None if compression == zipfile.ZIP_STORED else 1)
        archive.writestr(MANIFEST_NAME, manifest, compress_type=zipfile.ZIP_DEFLATED, compresslevel=1)
    checks = verify_zip(partial, rows, manifest)
    # Hashing source files again catches model/stage edits during compression,
    # even if the archive happened to contain a self-consistent older version.
    for item, record in zip(files, rows):
        if sha(item['path']) != record['sha256']:
            raise AssertionError('Source changed while packaging: ' + item['member'])
    partial.rename(destination)
    mp = destination.with_name(destination.stem + '_SHA256.tsv')
    mp.write_bytes(manifest)
    return {'filename': destination.name, 'bytes': destination.stat().st_size, 'sha256': sha(destination),
            'manifest_filename': mp.name, 'manifest_sha256': hashlib.sha256(manifest).hexdigest(),
            'uncompressed_payload_bytes': sum(x['bytes'] for x in rows), **checks}


def validate_stage(stage, report_dir):
    validator = stage / 'work/repro_validate.py'
    if not validator.is_file():
        raise FileNotFoundError('Stage the current validation entry point before packaging')
    # Require this reviewed validator rather than a stale staged copy.
    if sha(validator) != sha(WORK / 'repro_validate.py'):
        raise AssertionError('Staged validator differs from current reviewed source; rerun staging')
    report_dir.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, '-B', str(validator), '--require-complete', '--report-dir', str(report_dir)]
    with (report_dir / 'release_validation.log').open('w', encoding='utf-8') as stream:
        result = subprocess.run(command, cwd=stage, stdout=stream, stderr=subprocess.STDOUT, check=False)
    report = json.loads((report_dir / 'repro_validation_audit.json').read_text(encoding='utf-8'))
    if result.returncode or report['exit_code'] or report['n_failed'] or report['pending_stages'] or not report['require_complete']:
        raise AssertionError('The staged release failed complete validation; inspect ' + str(report_dir))
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--execute', action='store_true', help='Package after a fresh complete validation of the existing stage.')
    ap.add_argument('--stage-dir', type=Path, default=DEFAULT_STAGE)
    ap.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument('--with-corrected-matrices', action='store_true')
    args = ap.parse_args()
    stage, output = args.stage_dir.resolve(), args.output_dir.resolve()
    if not inside(stage, BASE) or not inside(output, BASE) or inside(output, stage):
        raise ValueError('Stage/output must be within this release workspace, with output outside the analysis stage')
    files = stage_files(stage)
    model, selected, pending = model_files(strict=False) if args.with_corrected_matrices else ([], [], [])
    plan = {'mode': 'execute' if args.execute else 'dry_run', 'stage': str(stage), 'stage_exists': stage.exists(),
            'output': str(output), 'stage_payload_file_count': len(files),
            'stage_payload_bytes': sum(x['path'].stat().st_size for x in files),
            'matrix_package_requested': args.with_corrected_matrices,
            'matrix_source_libraries_available': len(selected), 'pending_model_sources': pending,
            'matrix_payload_file_count_available': len(model),
            'stage_is_created_by_this_script': False, 'models_are_fitted_by_this_script': False,
            'execute_requires_fresh_complete_staged_validation': True,
            'archives': ['code_results.zip'] + (['corrected_matrices.zip'] if args.with_corrected_matrices else []),
            'compression': 'ZIP_STORED for already-compressed formats; deflate level 1 otherwise',
            'hash_scope': 'Every payload file in a SHA256 TSV; the manifest itself is hashed externally in release_integrity.json.',
            'self_reference_exclusions': [MANIFEST_NAME, INTEGRITY_NAME],
            'runtime_byproducts_excluded': ['__pycache__', '*.pyc', '*.pyo'],
            'duplicate_run_backups_excluded': ['work/repro_cached_runs'],
            'cached_reconstruction_evidence': 'outputs/reproducibility_v2/validation/',
            'publication_action': 'none'}
    print(json.dumps(plan, indent=2), flush=True)
    if not args.execute:
        return 0
    if not stage.exists() or not files:
        raise FileNotFoundError('The analysis stage does not exist; run the staging script after all upstream work completes')
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Choose a fresh release output directory: ' + str(output))
    output.mkdir(parents=True, exist_ok=True)
    validation_dir = stage / 'outputs/reproducibility_v2/validation/release_complete'
    validation = validate_stage(stage, validation_dir)
    # The new validation report/log are part of the finalized code/results payload.
    files = stage_files(stage)
    if args.with_corrected_matrices:
        model, selected, pending = model_files(strict=True)
        with (stage / 'outputs/reproducibility_v2/tables/ambient_selected_sources.tsv').open(encoding='utf-8', newline='') as stream:
            staged_sources = {(r['dataset'], r['sample']): r for r in csv.DictReader(stream, delimiter='\t')}
        for choice in selected:
            recorded = staged_sources[choice['dataset'], choice['sample']]
            for field in ['selected_audit_sha256', 'default_audit_sha256', 'selected_relative_directory']:
                if choice[field] != recorded[field]:
                    raise AssertionError('Matrix source changed since code/results staging: ' + choice['sample'] + '/' + field)
    packages = []
    packages.append(write_zip(output / 'code_results.zip', files, record_files(files)))
    if args.with_corrected_matrices:
        packages.append(write_zip(output / 'corrected_matrices.zip', model, record_files(model)))
    integrity = {'status': 'complete', 'created_utc': datetime.now(timezone.utc).isoformat(),
                 'packages': packages, 'complete_staged_validation': {
                     'n_checks': validation['n_checks'], 'n_failed': validation['n_failed'],
                     'pending_stages': validation['pending_stages'],
                     'audit_sha256': sha(validation_dir / 'repro_validation_audit.json'),
                     'validator_sha256': sha(stage / 'work/repro_validate.py')},
                 'selected_corrected_sources': selected, 'matrix_original_nonconverged_defaults_retained':
                     [x['sample'] for x in selected if x['numerical_extension_selected']],
                 'matrix_market_duplicate_not_included': True,
                 'manifest_does_not_hash_itself': True, 'integrity_record_is_external_to_archives': True,
                 'input_scope': ['explicit analysis stage', '11 named convergence-rule-selected model folders and their original default fits'],
                 'runtime_or_account_files_included': False, 'code_script_sha256': sha(Path(__file__)),
                 'public_deposition_performed': False}
    (output / INTEGRITY_NAME).write_text(json.dumps(integrity, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    print(json.dumps(integrity, indent=2), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
