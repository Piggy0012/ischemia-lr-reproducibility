"""Guarded v2 cached-table rebuild; defaults to a plan, never raw-cell inference.

Use --execute after inspecting the plan. --figures additionally creates PNG
previews; publication exports and visual approval are separate operations.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from repro_validate import BASE, WORK, TABLES, sha, validate

OUT = BASE / 'outputs/reproducibility_v2'
DEFAULT_STAGES = ['null', 'rank', 'raw_rank_summary', 'upstream', 'background', 'conditional', 'sorted', 'third_custom']
OPTIONAL_STAGES = ['ambient_summary', 'corrected_rank_summary', 'ambient_rank_matched', 'third_rank_summary']
DESCRIPTIONS = {
    'null': 'Rebuild whole-animal exact nulls from saved mean/fraction aggregates and original per-animal LIANA target scores.',
    'rank': 'Rebuild fixed-target percentile and rank-context diagnostics from saved scores and audit JSON.',
    'raw_rank_summary': 'Summarize saved raw full-network target exports; process()/load_sample() are forbidden.',
    'upstream': 'Run only the archived fixed-commit upstream aggregation functions on the unchanged raw full-network scores.',
    'background': 'Reaggregate saved full-network scores on a fixed global edge intersection; no upstream inference.',
    'conditional': 'Select the 20/10 conditional subsets of the existing 200 allocation statistics.',
    'sorted': 'Rebuild paired-pool descriptive effects/intervals from the author-normalized table; keep all old tests.',
    'third_custom': 'Rebuild third-cohort coavailability using six existing mean/fraction aggregates; streaming fallback forbidden.',
    'ambient_summary': 'Optional: summarize completed fractional-count expression effects and contamination tables; never fit DecontX.',
    'corrected_rank_summary': 'Optional: summarize 22 already-completed corrected-network targets; never run LIANA inference.',
    'ambient_rank_matched': 'Optional: compare raw/corrected scores on the same all-library/all-metric candidate intersection; no new P values.',
    'third_rank_summary': 'Optional: summarize six already-completed third-cohort networks and original comparison targets.',
    'core_figures': 'Optional: Figure 1/2 PNG previews from audited aggregate/source tables.',
    'rank_figure': 'Optional: Figure 3 PNG and grayscale previews, with a new pending-visual-review record.',
}
EXPECTED_REBUILD_VERSIONS = {'numpy': '2.3.5', 'pandas': '2.3.3', 'scipy': '1.18.1', 'liana': '1.10.0'}


def forbid(*args, **kwargs):
    raise RuntimeError('Raw-cell/model inference is forbidden by the cached-pipeline entry point')


def install_worker_guards(figures=False):
    """A missing cache cannot silently trigger a matrix read, download or fit."""
    def guard(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            value = str(args[0]).replace('\\', '/').lower()
            if value.endswith(('.npz', '.npy', '.h5ad', '.h5', '.mtx', '.mtx.gz', '.loom')):
                raise RuntimeError('Matrix access forbidden in cached worker: ' + value)
        if event in {'socket.connect', 'socket.getaddrinfo'}:
            raise RuntimeError('Network access forbidden in cached worker')
    sys.addaudithook(guard)
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
        os.environ[name] = '1'
    if not figures:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib.figure import Figure
        # The existing sorted script combines table and plot creation. Suppress
        # only figure writes, leaving its original statistical code untouched.
        Figure.savefig = lambda *args, **kwargs: None


def worker(stage, figures=False):
    install_worker_guards(figures)
    sys.argv = [str(Path(__file__))]  # Do not leak orchestration flags to imported scripts.
    if stage in {'null', 'rank', 'upstream', 'background', 'conditional', 'ambient_summary', 'ambient_rank_matched'}:
        names = {'null': 'repro_null_baseline', 'rank': 'repro_rank_mechanism',
                 'background': 'repro_rank_background', 'conditional': 'repro_conditional_null',
                 'upstream': 'repro_upstream_correction_audit',
                 'ambient_summary': 'repro_ambient_summary', 'ambient_rank_matched': 'repro_ambient_rank_comparison'}
        importlib.import_module(names[stage]).main()
    elif stage in {'raw_rank_summary', 'corrected_rank_summary'}:
        # Importing the module defines the inference function; it does not run it.
        m = importlib.import_module('repro_liana_rank_diagnostic')
        m.process = m.load_sample = forbid
        assert m.summarize('raw' if stage == 'raw_rank_summary' else 'decontx', ['primary', 'reference_singlet'])
    elif stage == 'sorted':
        runpy.run_path(str(WORK / 'repro_sorted_descriptive.py'), run_name='__main__')
    elif stage == 'third_custom':
        m = importlib.import_module('repro_third_coavailability')
        records = json.loads((m.ROOT / 'download_manifest.json').read_text(encoding='utf-8'))
        assert len(records) == 6
        audited = {}
        for rec in records:
            sample = rec['sample']
            for suffix in ['_stream_aggregation.json', '_mean_logcp10k.tsv.gz', '_fractions.tsv.gz', '_sample_counts.tsv']:
                p = m.FOLDER / (sample + suffix)
                if not p.is_file():
                    raise FileNotFoundError('Required cached aggregate missing; no streaming fallback: ' + str(p))
            a = json.loads((m.FOLDER / (sample + '_stream_aggregation.json')).read_text())
            assert a['matrix_gene_major_unique'] and a['all_assigned_counts_conserved']
            audited[sample] = a
        m.aggregate = lambda rec: audited[rec['sample']]
        m.chunks = forbid
        # Metadata amendments are retained separately from the recomputed statistics.
        ap = TABLES / 'third_coavailability_audit.json'
        prior = json.loads(ap.read_text())
        m.main()
        new = json.loads(ap.read_text())
        for key in ['author_reported_biological_replicates_per_condition', 'biological_replicate_source']:
            if key in prior:
                new[key] = prior[key]
        ap.write_text(json.dumps(new, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    elif stage == 'third_rank_summary':
        m = importlib.import_module('repro_third_liana')
        m.process = m.build_cache = m.load_context = forbid
        m.summarize()
    elif stage == 'core_figures':
        m = importlib.import_module('repro_core_figures')
        m.EXPORT = False
        m.figure1(); m.figure2()
    elif stage == 'rank_figure':
        m = importlib.import_module('repro_rank_figure')
        m.EXPORT = False
        m.build()
    else:
        raise ValueError(stage)


def output_paths():
    return sorted(p for folder in ['tables', 'figure_source_data', 'figures']
                  for p in (OUT / folder).rglob('*') if p.is_file())


def protected_paths():
    """Historical tests and biological inputs are never outputs of this pipeline."""
    files = set()
    files.add(WORK / 'literature/mouseconsensus.csv')
    for folder in ['results', 'revision_analysis/results', 'sorted_rna', 'processed', 'revision_analysis/processed',
                   'repro_liana_rank_diagnostic/raw', 'repro_liana_rank_diagnostic/decontx',
                   'repro_decontx_analysis/processed', 'repro_decontx_analysis/results',
                   'repro_third_liana', 'repro_third_cohort/processed']:
        for p in (WORK / folder).rglob('*'):
            if not p.is_file() or not p.name.endswith(('.tsv', '.tsv.gz', '.json')):
                continue
            # The third script rewrites this derived concatenation; preserve it in
            # the output backup and compare semantically instead of byte hashing.
            if p == WORK / 'repro_third_cohort/processed/GSE332910/sample_cell_counts.tsv':
                continue
            # Per-cell annotation matrices/large model outputs are not consumed.
            if p.name.endswith(('_cells.tsv.gz', '_features.tsv.gz')):
                continue
            files.add(p)
    return sorted(files)


def semantic_compare(before, after):
    if not after.exists():
        return {'status': 'missing_after', 'scientific_values_equal': False}
    if not before.exists():
        return {'status': 'new_output', 'scientific_values_equal': None}
    first, second = sha(before), sha(after)
    result = {'before_sha256': first, 'after_sha256': second}
    if after.name.endswith(('.tsv', '.tsv.gz', '.csv')):
        sep = ',' if after.suffix == '.csv' else '\t'
        a, b = pd.read_csv(before, sep=sep, float_precision='round_trip'), pd.read_csv(after, sep=sep, float_precision='round_trip')
        result.update(before_rows=len(a), after_rows=len(b))
        if list(a.columns) != list(b.columns) or a.shape != b.shape:
            return result | {'status': 'schema_or_row_count_changed', 'scientific_values_equal': False}
        changed, maxdiff = [], 0.0
        for col in a:
            av, bv = a[col], b[col]
            if not np.array_equal(av.isna().to_numpy(), bv.isna().to_numpy()):
                changed.append(col); continue
            if pd.api.types.is_numeric_dtype(av) and pd.api.types.is_numeric_dtype(bv):
                finite = av.notna()
                if finite.any():
                    x, y = av[finite].to_numpy(), bv[finite].to_numpy()
                    diff = float(np.max(np.abs(x.astype(float) - y.astype(float))))
                    maxdiff = max(maxdiff, diff)
                    integer = pd.api.types.is_integer_dtype(av) and pd.api.types.is_integer_dtype(bv)
                    equal = np.array_equal(x, y) if integer else np.allclose(x, y, atol=1e-12, rtol=1e-10)
                    if not equal:
                        changed.append(col)
            elif not av.fillna('<MISSING>').astype(str).equals(bv.fillna('<MISSING>').astype(str)):
                changed.append(col)
        return result | {'status': 'semantic_equal' if not changed else 'table_values_changed',
                         'scientific_values_equal': not changed, 'changed_columns': '|'.join(changed),
                         'max_numeric_absolute_difference': maxdiff}
    if after.suffix == '.json':
        a, b = json.loads(before.read_text(encoding='utf-8')), json.loads(after.read_text(encoding='utf-8'))
        return result | {'status': 'json_equal' if a == b else 'audit_metadata_changed',
                         'scientific_values_equal': None,
                         'note': 'Source hashes, paths and review status are retained; inspect changed metadata and post-validation.'}
    return result | {'status': 'bytes_equal' if first == second else 'binary_or_text_changed',
                     'scientific_values_equal': None}


def build_plan(args, validation):
    stages = list(args.stages or DEFAULT_STAGES)
    skipped = []
    pending = {r['stage'] for r in validation['pending_stages']}
    if args.include_completed_optional:
        for stage in OPTIONAL_STAGES:
            needed = 'third_liana' if stage == 'third_rank_summary' else 'ambient'
            if needed in pending:
                skipped.append({'stage': stage, 'reason': 'Completed upstream audit absent; no model is started.'})
            else:
                stages.append(stage)
    if args.figures:
        stages += ['core_figures', 'rank_figure']
    return {'stages': [{'stage': x, 'operation': DESCRIPTIONS[x]} for x in stages], 'skipped': skipped,
            'execute_requested': args.execute, 'figures_requested': args.figures,
            'python': sys.executable, 'working_directory': str(WORK),
            'expected_rebuild_versions': EXPECTED_REBUILD_VERSIONS,
            'model_fit_stages': [], 'network_access': False, 'raw_matrix_access': False,
            'publication_figure_exports': False,
            'visual_review': 'Any regenerated PNG must be visually inspected; earlier PDF/SVG exports are retained, not silently reapproved.',
            'pending_stages': validation['pending_stages']}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--execute', action='store_true', help='Run the guarded cache-only stages; otherwise print the plan.')
    ap.add_argument('--stages', nargs='+', choices=DEFAULT_STAGES)
    ap.add_argument('--include-completed-optional', action='store_true')
    ap.add_argument('--figures', action='store_true', help='Also regenerate PNG previews; never auto-export/reapprove PDF/SVG.')
    ap.add_argument('--require-complete', action='store_true', help='Require final ambient and third-LIANA completion audits.')
    ap.add_argument('--run-dir', type=Path)
    ap.add_argument('--worker', choices=list(DESCRIPTIONS), help=argparse.SUPPRESS)
    ap.add_argument('--worker-figures', action='store_true', help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.worker:
        worker(args.worker, args.worker_figures)
        return 0
    pre = validate(require_complete=args.require_complete)
    plan = build_plan(args, pre)
    print(json.dumps(plan, indent=2), flush=True)
    if pre['exit_code']:
        print('Preflight validation failed or required stages are incomplete.', flush=True)
        for c in pre['checks']:
            if not c['passed']:
                print('FAIL', c['check'], c['detail'])
        return pre['exit_code']
    if not args.execute:
        return 0
    mismatched = {x: {'expected': wanted, 'installed': importlib.metadata.version(x)}
                  for x, wanted in EXPECTED_REBUILD_VERSIONS.items() if importlib.metadata.version(x) != wanted}
    if mismatched:
        print('Rebuild environment differs from the frozen analysis environment: ' + json.dumps(mismatched), flush=True)
        return 1
    run = args.run_dir or WORK / 'repro_cached_runs' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    if run.exists():
        raise FileExistsError('Choose a new run directory so previous evidence is never overwritten: ' + str(run))
    run.mkdir(parents=True)
    before = run / 'before'
    before.mkdir()
    (run / 'plan.json').write_text(json.dumps(plan, indent=2), encoding='utf-8')
    (run / 'validation_before.json').write_text(json.dumps(pre, indent=2, ensure_ascii=False), encoding='utf-8')
    original_outputs = output_paths()
    derived_count = WORK / 'repro_third_cohort/processed/GSE332910/sample_cell_counts.tsv'
    original_outputs += [derived_count]
    for path in original_outputs:
        dest = before / path.relative_to(BASE)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    protected = {str(p.relative_to(BASE)): sha(p) for p in protected_paths()}
    (run / 'protected_inputs_sha256.json').write_text(json.dumps(protected, indent=2), encoding='utf-8')
    runs, error = [], None
    for item in plan['stages']:
        stage = item['stage']
        command = [sys.executable, str(Path(__file__).resolve()), '--worker', stage]
        if args.figures:
            command.append('--worker-figures')
        log = run / (stage + '.log')
        print('Cached stage:', stage, flush=True)
        started = time.monotonic()
        with log.open('w', encoding='utf-8') as stream:
            result = subprocess.run(command, cwd=WORK, stdout=stream, stderr=subprocess.STDOUT, check=False,
                                    env=os.environ | {'PYTHONUTF8': '1'})
        runs.append({'stage': stage, 'exit_code': result.returncode, 'elapsed_seconds': time.monotonic() - started,
                     'log': log.name, 'command_arguments': command})
        if result.returncode:
            error = 'Stage failed: ' + stage
            break
    comparisons = []
    after_paths = set(output_paths() + [derived_count])
    for path in sorted(set(original_outputs) | after_paths):
        old = before / path.relative_to(BASE)
        comparisons.append({'path': str(path.relative_to(BASE)), **semantic_compare(old, path)})
    changed_protected = [p for p, expected in protected.items() if not (BASE / p).exists() or sha(BASE / p) != expected]
    post = validate(require_complete=args.require_complete, report_dir=run / 'validation_after')
    table_changes = [x for x in comparisons if x.get('scientific_values_equal') is False]
    code = 1 if error or changed_protected or table_changes or post['exit_code'] == 1 else post['exit_code']
    report = {'status': 'passed' if code == 0 else 'failed_or_incomplete', 'exit_code': code,
              'stages': runs, 'error': error, 'protected_input_changes': changed_protected,
              'scientific_table_changes': table_changes, 'table_comparison_atol': 1e-12,
              'table_comparison_rtol': 1e-10, 'table_row_and_column_order_compared': True,
              'gzip_container_timestamps_ignored_for_semantics': True,
              'new_model_fits': False, 'matrix_reads_guarded': True, 'network_reads_guarded': True,
              'visual_review_required_for_regenerated_png': args.figures,
              'old_tests_preserved': not changed_protected, 'before_outputs_retained': 'before/',
              'post_validation_status': post['status'], 'pending_stages': post['pending_stages'],
              'versions': {x: importlib.metadata.version(x) for x in ['numpy', 'pandas', 'scipy', 'matplotlib', 'liana']},
              'script_sha256': sha(Path(__file__))}
    pd.DataFrame(comparisons).to_csv(run / 'semantic_comparison.tsv', sep='\t', index=False)
    (run / 'cached_pipeline_audit.json').write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    print(json.dumps({'status': report['status'], 'exit_code': code, 'run_directory': str(run),
                      'protected_input_changes': changed_protected, 'scientific_table_change_count': len(table_changes),
                      'pending_stages': post['pending_stages']}, indent=2), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
