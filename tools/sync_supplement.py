"""Copy an explicit reviewed supplement manifest. Default is dry run; never publish."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
BLOCKED = {'.git', '.ssh', '.aws', 'r_runtime', 'revision_env', 'venv', '.venv', 'rlibrary',
           'node_modules', 'site-packages', 'renv_tooling', 'renv_bootstrap', 'r-library', 'library',
           'runtime', 'r-runtime', 'repro_v3_github_cli', 'repro_cached_runs'}


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('manifest', type=Path)
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    records = json.loads(args.manifest.read_text(encoding='utf-8'))['files']
    prepared = []
    for row in records:
        relative = PurePosixPath(row['destination'])
        assert not relative.is_absolute() and '..' not in relative.parts and ':' not in str(relative) and '\\' not in str(relative)
        assert not set(relative.parts) & BLOCKED
        assert relative.suffix.lower() not in {'.zip', '.exe', '.dll', '.pyd', '.key', '.pem'}
        assert relative.name.lower() not in {'.env', '.netrc', 'auth.json', 'token.json', 'tokens.json', 'credentials.json', 'cookies.txt'}
        source, target = Path(row['source']).resolve(), ROOT.joinpath(*relative.parts)
        assert source.is_file() and not source.is_symlink() and sha(source) == row['sha256']
        assert source.stat().st_size < 100 * 2**20
        if target.exists() and sha(target) != row['sha256']:
            assert row.get('expected_existing_sha256') == sha(target), 'Replacement requires the exact reviewed old hash.'
        prepared.append((source, target, row))
    print(json.dumps({'mode': 'execute' if args.execute else 'dry_run', 'files': [r for _, _, r in prepared],
                      'external_publication': False}, indent=2))
    if args.execute:
        for source, target, row in prepared:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            assert sha(target) == row['sha256']
        subprocess.run([sys.executable, str(ROOT / 'tools/update_file_index.py')], check=True, cwd=ROOT)


if __name__ == '__main__':
    main()
