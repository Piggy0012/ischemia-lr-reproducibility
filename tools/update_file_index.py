"""Index publication bytes and verify the figure/source closure; no Git operations."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parent.parent
EXCLUDE = {'.git', '__pycache__', 'repro_cached_runs', 'revision_env', 'venv', '.venv', 'r_runtime',
           'rlibrary', 'site-packages', 'renv_tooling', 'renv_bootstrap', 'r-library', 'node_modules'}
SELF_FILES = {'FILE_INDEX.tsv', 'FILE_INDEX_SHA256.txt'}
FORBIDDEN_NAMES = {'.env', '.netrc', '.npmrc', 'credentials.json', 'credentials.toml', 'auth.json',
                   'token.json', 'tokens.json', 'cookies.txt', 'id_rsa', 'id_ed25519'}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def eligible():
    result = []
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or any(x in EXCLUDE or x.startswith('restore_smoke_library') for x in p.relative_to(ROOT).parts):
            continue
        if p.relative_to(ROOT).as_posix() in SELF_FILES or p.suffix in {'.pyc', '.pyo'}:
            continue
        assert not p.is_symlink(), 'External/symbolic source is not allowed: ' + str(p)
        assert p.name.lower() not in FORBIDDEN_NAMES and not p.name.startswith('.env.'), str(p)
        assert p.suffix.lower() not in {'.zip', '.exe', '.dll', '.pyd', '.msi', '.pem', '.key', '.npz'}, str(p)
        result.append(p)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    files = eligible()
    oversized = [p.relative_to(ROOT).as_posix() for p in files if p.stat().st_size >= 100 * 2**20]
    assert not oversized, oversized
    with (ROOT / 'FIGURE_SOURCE_INDEX.tsv').open(encoding='utf-8', newline='') as f:
        figure_rows = list(csv.DictReader(f, delimiter='\t'))
    assert len({r['figure_id'] for r in figure_rows}) == 10
    for row in figure_rows:
        p = ROOT / row['path']
        assert p.is_file() and p.stat().st_size == int(row['bytes']) and sha(p) == row['sha256'], row['path']
    if args.verify:
        with (ROOT / 'FILE_INDEX.tsv').open(encoding='utf-8', newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        assert {r['path'] for r in rows} == {p.relative_to(ROOT).as_posix() for p in files}
        for row in rows:
            p = ROOT / row['path']
            assert p.stat().st_size == int(row['bytes']) and sha(p) == row['sha256'], row['path']
        expected = (ROOT / 'FILE_INDEX_SHA256.txt').read_text().split()[0]
        assert sha(ROOT / 'FILE_INDEX.tsv') == expected
    else:
        rows = [{'path': p.relative_to(ROOT).as_posix(), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in files]
        with (ROOT / 'FILE_INDEX.tsv').open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, ['path', 'bytes', 'sha256'], delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows)
        (ROOT / 'FILE_INDEX_SHA256.txt').write_text(sha(ROOT / 'FILE_INDEX.tsv') + '  FILE_INDEX.tsv\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'mode': 'verify' if args.verify else 'write', 'indexed_files': len(rows),
                      'payload_bytes': sum(int(r['bytes']) for r in rows), 'figure_count': 10,
                      'figure_source_rows': len(figure_rows), 'files_at_least_100_MiB': oversized,
                      'index_self_reference_exclusions': sorted(SELF_FILES)}, indent=2))


if __name__ == '__main__':
    main()
