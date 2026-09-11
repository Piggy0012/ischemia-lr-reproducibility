"""Capture the actual environment and release sizes; never mutate analysis outputs."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import importlib
import importlib.metadata as md
import json
import platform
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
WORK = HERE.parent
BASE = WORK.parent


def main():
    assert sys.version_info[:3] == (3, 12, 14)
    pip_check = subprocess.run([sys.executable, '-m', 'pip', 'check'], capture_output=True, text=True)
    frozen = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze', '--all'], text=True)
    # The host-only artifact package is neither imported nor required by any
    # scientific reconstruction. Its machine-local URL is intentionally excluded.
    retained, excluded = [], []
    for line in frozen.splitlines():
        if line.lower().startswith(('artifact-tool-v2', 'artifact_tool_v2')):
            excluded.append({'package': 'artifact-tool-v2', 'reason': 'Host artifact tooling; local path dependency; unused by the scientific scripts.'})
        elif line.strip():
            assert '==' in line and ' @ ' not in line, 'Unresolved non-version dependency: ' + line
            retained.append(line)
    retained.sort(key=str.casefold)
    (HERE / 'requirements-python.lock.txt').write_text('\n'.join(retained) + '\n', encoding='utf-8')
    lines = ['# Observed Windows analysis environment; conda solver/restore not tested.',
             '# Pip package versions are exact. Conda build strings are not locked.',
             'name: ischemia-lr-repro', 'channels:', '  - conda-forge', '  - nodefaults',
             'dependencies:', '  - python=3.12.14', '  - pip=25.0.1', '  - pip:']
    lines += ['      - ' + x for x in retained if not x.lower().startswith('pip==')]
    (HERE / 'environment.yml').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    imports = {}
    for name in ['numpy', 'pandas', 'scipy', 'anndata', 'scanpy', 'liana', 'matplotlib', 'celltypist', 'pydeseq2']:
        module = importlib.import_module(name)
        imports[name] = {'metadata_version': md.version(name), 'imported_version': getattr(module, '__version__', None),
                         'module_file': str(Path(module.__file__).resolve())}
    (HERE / 'python_environment.json').write_text(json.dumps({
        'captured_utc': datetime.now(timezone.utc).isoformat(), 'python': platform.python_version(),
        'python_executable': sys.executable, 'platform': platform.platform(),
        'imports': imports, 'pip_check_exit_code': pip_check.returncode,
        'scrublet_implementation': 'scanpy.pp.scrublet (no separately installed scrublet distribution)',
        'pip_check_stdout': pip_check.stdout, 'pip_check_stderr': pip_check.stderr,
        'pip_package_count': len(retained), 'excluded': excluded,
        'conda_available': shutil.which('conda'), 'mamba_available': shutil.which('mamba'),
        'micromamba_available': shutil.which('micromamba'), 'fresh_conda_restore_tested': False,
        'current_environment_changed': False, 'new_model_fit': False}, indent=2), encoding='utf-8')
    rows = []
    release = WORK / 'repro_delivery/release'
    for p in sorted(release.glob('*')):
        if not p.is_file():
            continue
        n = p.stat().st_size
        rows.append({'scope': 'release_file', 'archive': '', 'path': p.relative_to(BASE).as_posix(),
                     'bytes': n, 'MiB': n / 2**20, 'over_100_MB': n > 100_000_000,
                     'over_100_MiB': n > 100 * 2**20, 'at_least_2_GiB': n >= 2 * 2**30})
        if p.suffix == '.zip':
            with zipfile.ZipFile(p) as z:
                for item in z.infolist():
                    if item.is_dir():
                        continue
                    n = item.file_size
                    rows.append({'scope': 'zip_member', 'archive': p.name, 'path': item.filename,
                                 'bytes': n, 'MiB': n / 2**20, 'over_100_MB': n > 100_000_000,
                                 'over_100_MiB': n > 100 * 2**20, 'at_least_2_GiB': n >= 2 * 2**30})
    for filename, selected in [('release_file_sizes.tsv', rows),
                                ('release_files_over_100MB.tsv', [r for r in rows if r['over_100_MB']])]:
        with (HERE / filename).open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(selected)
    print(json.dumps({'pip_check': pip_check.returncode, 'python_locked': len(retained),
                      'excluded_host_packages': len(excluded), 'size_records': len(rows),
                      'large_release_files': [r for r in rows if r['scope'] == 'release_file' and r['over_100_MB']],
                      'large_zip_members': sum(r['scope'] == 'zip_member' and r['over_100_MiB'] for r in rows)}, indent=2))


if __name__ == '__main__':
    main()
