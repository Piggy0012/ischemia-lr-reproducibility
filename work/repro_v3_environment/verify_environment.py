"""Verify version pins and the two distinct LIANA source scopes; no analysis writes."""
from pathlib import Path
import argparse
import ast
import hashlib
import importlib.metadata as md
import importlib.util
import json
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MERGE = 'd4211373692e7b9c10210488ccb1efe06452b097'
HASHES = {
    'liana_merge_d421137_aggregate.py': '484b2e61d7bb01d141116a6fdf8f8b7d269973396be3354d10702222770aa875',
    'liana_1_10_0_aggregate.py': '7e7fa25f940e1862db2082ad79a81ccadfa64ad317fe9fb05d928b52c9c7984b',
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--release-root', type=Path, default=HERE.parent.parent)
    ap.add_argument('--output', type=Path, default=HERE / 'environment_verification.json')
    args = ap.parse_args()
    work = args.release_root.resolve() / 'work'
    checks = []
    def check(name, passed, detail=None):
        checks.append({'check': name, 'passed': bool(passed), 'detail': detail})
    check('python_3.12.14', platform.python_version() == '3.12.14')
    for line in (HERE / 'requirements-python.lock.txt').read_text().splitlines():
        package, expected = line.split('==')
        try:
            actual = md.version(package)
        except md.PackageNotFoundError:
            actual = None
        check('pip_pin/' + package, actual == expected, {'expected': expected, 'actual': actual})
    for name, expected in HASHES.items():
        p = work / 'repro_literature' / name
        check('archived_source/' + name, p.is_file() and sha(p) == expected)
    spec = importlib.util.find_spec('liana')
    current = Path(spec.origin).parent / 'method/_pipe_utils/_aggregate.py'
    check('installed_1.10.0_aggregate_source', sha(current) == HASHES['liana_1_10_0_aggregate.py'])
    original = ast.parse((work / 'repro_literature/liana_merge_d421137_aggregate.py').read_text())
    extracted = ast.parse((work / 'repro_literature/liana_merge_d421137_extracted_aggregation.py').read_text())
    functions = {n.name: ast.dump(n, include_attributes=False) for n in original.body if isinstance(n, ast.FunctionDef)}
    included = {n.name: ast.dump(n, include_attributes=False) for n in extracted.body if isinstance(n, ast.FunctionDef)}
    for name, body in included.items():
        if name == '_logg':
            continue  # Deliberate logging-only stub, not an upstream aggregation function.
        check('fixed_commit_function_AST/' + name, functions.get(name) == body)
    audit_path = args.release_root.resolve() / 'outputs/reproducibility_v2/tables/rank_upstream_audit.json'
    previous = json.loads(audit_path.read_text())
    check('archived_22_network_scope', previous['n_library_config_audits'] == 22 and
          previous['upstream_merge_commit'] == MERGE and previous['all_upstream_vs_original_unique_diagnostic_bitwise_equal'] and
          previous['maximum_absolute_difference_from_original_unique_diagnostic'] == 0 and
          not previous['new_full_pipeline_run'])
    pip = subprocess.run([sys.executable, '-m', 'pip', 'check'], capture_output=True, text=True)
    check('pip_check', pip.returncode == 0, pip.stdout + pip.stderr)
    report = {'status': 'passed' if all(r['passed'] for r in checks) else 'failed', 'checks': checks,
              'n_checks': len(checks), 'new_full_network_fit': False, 'new_22_network_aggregation_replay': False,
              'existing_22_network_audit_sha256': sha(audit_path), 'fixed_commit': MERGE,
              'fresh_environment_restore_tested': False, 'python_executable': sys.executable,
              'script_sha256': sha(Path(__file__))}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: report[k] for k in ['status', 'n_checks', 'new_full_network_fit', 'new_22_network_aggregation_replay']}, indent=2))
    return int(report['status'] != 'passed')


if __name__ == '__main__':
    sys.exit(main())
