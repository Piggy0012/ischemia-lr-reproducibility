"""Preserve native LIANA full networks and diagnose duplicate-column rank behavior.

CLI: python work/repro_liana_rank_diagnostic.py --count-source raw --configs primary reference_singlet
Callable: process(acc, sample, config='primary', count_source='raw').
The original package is never modified. Diagnostic ranks are named explicitly.
"""
from __future__ import annotations
import argparse
import gc
import hashlib
import importlib.metadata
import io
import json
import os
import time
from pathlib import Path
for _n in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ[_n] = '1'
import anndata as ad
import liana as li
import numpy as np
import pandas as pd
from scipy import sparse, stats
from liana.method._pipe_utils._aggregate import _rank_aggregate
from revision_data import load_sample, samples
from repro_decontx_source import resolve

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'repro_liana_rank_diagnostic'
TABLES = ROOT.parent / 'outputs/reproducibility_v2/tables'
KEY = ['source', 'target', 'ligand_complex', 'receptor_complex']
RENAMES = dict(zip(KEY, ['sender', 'receiver', 'ligand', 'receptor']))
DATASETS = ['GSE174574', 'GSE245386']
METRICS = ['lr_means', 'expr_prod', 'lrscore', 'native_magnitude_priority', 'diagnostic_unique_column_priority']
TOL = 1e-12

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def stems(acc, sample, config, count_source):
    return OUT / count_source / acc / (config + '__' + sample)

def process(acc, sample, config='primary', count_source='raw', force=False):
    stem = stems(acc, sample, config, count_source)
    audit_path = Path(str(stem) + '.json')
    full_path = Path(str(stem) + '__full_network.tsv.gz')
    target_path = Path(str(stem) + '__target.tsv.gz')
    source_selection = None
    if count_source == 'decontx':
        corrected_folder, source_selection = resolve(acc, sample)
    if not force and all(p.exists() for p in [audit_path, full_path, target_path]):
        audit = json.loads(audit_path.read_text(encoding='utf-8'))
        assert audit['status'] == 'complete' and audit['count_source'] == count_source
        print('DIAGNOSTIC EXISTS', count_source, acc, sample, config, flush=True)
        if count_source == 'raw' or audit.get('corrected_source_selection') == source_selection:
            return audit
    started = time.monotonic()
    x, q, genes = load_sample(acc, sample)
    corrected_hash = None
    corrected_order_verified = None
    if count_source == 'decontx':
        corrected_path = corrected_folder / 'corrected_counts.npz'
        corrected_hash = sha(corrected_path)
        corrected_audit = json.loads((corrected_path.parent / 'audit.json').read_text(encoding='utf-8'))
        assert corrected_audit['output_sha256']['corrected_counts.npz'] == corrected_hash
        corrected_cells = pd.read_csv(corrected_path.parent / 'cells.tsv.gz', sep='\t', usecols=['barcode'])
        corrected_genes = pd.read_csv(corrected_path.parent / 'genes.tsv.gz', sep='\t', usecols=['symbol'])
        assert np.array_equal(corrected_cells.barcode.astype(str), q.barcode.astype(str)), 'Corrected barcode order differs'
        assert np.array_equal(corrected_genes.symbol.astype(str), genes.astype(str)), 'Corrected gene order differs'
        corrected_order_verified = True
        corrected = sparse.load_npz(corrected_path).tocsr()
        assert corrected.shape == x.shape, 'Corrected CSR must retain every raw QC cell and gene in original order'
        assert np.isfinite(corrected.data).all() and np.all(corrected.data >= 0)
        del x; x = corrected
    else:
        assert count_source == 'raw'
    totals = np.asarray(x.sum(axis=1)).ravel().astype(float)
    assert np.isfinite(totals).all() and np.all(totals > 0), 'Zero/nonfinite corrected cell totals require explicit handling'
    if count_source == 'raw':
        assert np.array_equal(totals, q.n_umis.to_numpy()), 'Raw all-gene totals do not reproduce original totals'
    keep = (q.cell_type != 'Unassigned').to_numpy()
    identity_hash = None
    if config == 'reference_singlet':
        identity_path = ROOT / f'revision_results/identity/{acc}/{sample}_identity.tsv.gz'
        r = pd.read_csv(identity_path, sep='\t')
        identity_hash = sha(identity_path)
        assert np.array_equal(q.barcode, r.barcode)
        target = q.cell_type.isin(['Astrocyte', 'Endothelial']).to_numpy()
        supported = (r.whole_brain_broad == q.cell_type).to_numpy() & (r.whole_brain_probability.to_numpy() >= .5)
        keep &= ~r.predicted_doublet.to_numpy(bool) & (~target | supported)
    else:
        assert config == 'primary'
    group_counts = q.loc[keep, 'cell_type'].value_counts()
    retained = group_counts[group_counts >= 30].index
    keep &= q.cell_type.isin(retained).to_numpy()
    assert all(t in retained for t in ['Astrocyte', 'Endothelial'])
    sub = q.loc[keep]
    y = x[keep]
    n_positive_before_float32 = int(np.count_nonzero(y.data > 0))
    y = y.astype(np.float32)
    del x; gc.collect()
    y = y.multiply((1e4 / totals[keep]).astype(np.float32)[:, None]).tocsr()
    y.data = np.log1p(y.data)
    # LIANA 1.10.0 _get_props counts stored sparse entries with getnnz.
    # Remove explicit zeros after normalisation so underflowed fractional
    # estimates cannot be counted as detected expression.
    n_explicit_zeros_removed = int(np.count_nonzero(y.data == 0))
    n_positive_after_float32 = int(np.count_nonzero(y.data > 0))
    y.eliminate_zeros()
    a = ad.AnnData(y, obs=pd.DataFrame({'cell_type': pd.Categorical(sub.cell_type.to_numpy())},
                   index=sub.barcode.astype(str)), var=pd.DataFrame(index=genes))
    resource_path = ROOT / 'literature/mouseconsensus.csv'
    resource = pd.read_csv(resource_path)[['source_genesymbol', 'target_genesymbol']].drop_duplicates().rename(
        columns={'source_genesymbol': 'ligand', 'target_genesymbol': 'receptor'})
    print('DIAGNOSTIC START', count_source, acc, sample, config, a.shape, flush=True)
    result = li.mt.rank_aggregate(a, groupby='cell_type', resource=resource, expr_prop=.1, min_cells=30,
        use_raw=False, n_perms=100, seed=20260910, n_jobs=1, return_all_lrs=False, inplace=False, verbose=False)
    assert not result.duplicated(KEY).any()
    native_specs = li.mt.rank_aggregate.magnitude_specs
    unique_specs, used = {}, set()
    for name, spec in native_specs.items():
        if spec[0] not in used:
            unique_specs[name] = spec; used.add(spec[0])
    # Both calculations use the same full-network rows, never a target-only rank universe.
    native_recomputed = _rank_aggregate(result.copy(), native_specs, 'rra')
    native_error = float(np.max(np.abs(native_recomputed - result.magnitude_rank.to_numpy())))
    assert native_error < TOL, f'Cannot reproduce native rank from complete full network: {native_error}'
    result['diagnostic_unique_column_magnitude_rank'] = _rank_aggregate(result.copy(), unique_specs, 'rra')
    result['native_magnitude_priority'] = 1 - result.magnitude_rank
    result['diagnostic_unique_column_priority'] = 1 - result.diagnostic_unique_column_magnitude_rank
    target_mask = (((result.source == 'Astrocyte') & (result.target == 'Endothelial')) |
                   ((result.source == 'Endothelial') & (result.target == 'Astrocyte')))
    target_result = result.loc[target_mask].copy()
    condition = 'Sham' if ('sham' in q.prefix.iloc[0].lower() or '_WTC' in q.prefix.iloc[0]) else 'MCAO'
    prior_errors = {}
    old_ap = ROOT / f'revision_results/liana/{acc}/{config}__{sample}.json'
    old_tp = Path(str(old_ap).replace('.json', '.tsv.gz'))
    if count_source == 'raw':
        old = pd.read_csv(old_tp, sep='\t').set_index(KEY).sort_index()
        # The archived TSV formats float32 scores as their shortest decimal
        # representations. Compare after the identical TSV round-trip rather
        # than comparing parsed decimals to binary float32 promoted to float64.
        now = pd.read_csv(io.StringIO(target_result.to_csv(sep='\t', index=False)), sep='\t').set_index(KEY).sort_index()
        assert old.index.equals(now.index), 'Raw rerun eligible target keys differ from saved original'
        numeric = old.select_dtypes(include=[np.number]).columns
        for col in numeric:
            error = float(np.nanmax(np.abs(old[col].to_numpy() - now[col].to_numpy())))
            assert error < 1e-9, f'Raw rerun mismatch in {col}: {error}'
            prior_errors[col] = error
        old_audit = json.loads(old_ap.read_text(encoding='utf-8'))
        assert old_audit['all_pair_rows'] == len(result)
        assert old_audit['context_counts'] == sub.cell_type.value_counts().to_dict()
    for df in [result, target_result]:
        df.insert(0, 'sample', sample); df.insert(0, 'dataset', acc)
        df.insert(0, 'config', config); df.insert(0, 'count_source', count_source)
        df['condition'] = condition; df['expression_eligible'] = True
    stem.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(full_path, sep='\t', index=False)
    target_result.to_csv(target_path, sep='\t', index=False)
    import liana.method._pipe_utils._aggregate as aggregate_module
    report = {'status': 'complete', 'count_source': count_source, 'dataset': acc, 'sample': sample,
        'corrected_source_selection': source_selection,
        'normalised_matrix_dtype': str(y.dtype),
        'n_positive_before_float32': n_positive_before_float32,
        'n_positive_after_float32': n_positive_after_float32,
        'n_explicit_zeros_removed': n_explicit_zeros_removed,
        'sparse_detection_audit': 'Explicit zeros removed before LIANA getnnz-based expression proportions',
        'config': config, 'condition': condition, 'liana_version': importlib.metadata.version('liana'),
        'methods': [m.method_name for m in li.mt.rank_aggregate.methods],
        'native_magnitude_specs': native_specs, 'diagnostic_unique_score_specs': unique_specs,
        'parameters': {'expr_prop': .1, 'min_cells': 30, 'use_raw': False, 'n_perms': 100, 'seed': 20260910,
                       'n_jobs': 1, 'return_all_lrs': False, 'aggregate_method': 'rra'},
        'computed_all_cell_type_pairs': True, 'all_pair_rows': len(result), 'exported_target_rows': len(target_result),
        'n_input_cells': len(q), 'n_context_cells': a.n_obs, 'n_context_genes': a.n_vars,
        'context_counts': sub.cell_type.value_counts().to_dict(),
        'denominator': 'per-cell total of supplied count source across all genes, before cell-type filtering',
        'zero_total_cells': int((totals == 0).sum()), 'corrected_counts_sha256': corrected_hash,
        'corrected_barcode_gene_order_verified': corrected_order_verified,
        'identity_sha256': identity_hash, 'resource_sha256': sha(resource_path),
        'barcode_order_sha256': hashlib.sha256('\n'.join(q.barcode.astype(str)).encode()).hexdigest(),
        'gene_order_sha256': hashlib.sha256('\n'.join(genes.astype(str)).encode()).hexdigest(),
        'source_aggregate_sha256': sha(aggregate_module.__file__), 'script_sha256': sha(__file__),
        'plan_sha256': sha(ROOT / 'repro_liana_rank_diagnostic_plan.md'),
        'native_full_network_reconstruction_max_abs_error': native_error,
        'raw_target_scores_verified_against_prior': count_source == 'raw',
        'raw_prior_comparison_representation': 'both scores parsed from pandas TSV decimal serialization',
        'raw_target_prior_max_abs_errors': prior_errors,
        'prior_target_sha256': sha(old_tp) if count_source == 'raw' else None,
        'full_network_sha256': sha(full_path), 'target_sha256': sha(target_path),
        'package_modified': False, 'diagnostic_is_official_package_release': False,
        'elapsed_seconds': time.monotonic() - started}
    audit_path.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print('DIAGNOSTIC DONE', count_source, sample, config, 'network', len(result),
          'seconds', round(report['elapsed_seconds'], 1), flush=True)
    del result, target_result, a, y; gc.collect()
    return report

def summarize(count_source, configs):
    all_rows, effects, means, tie_rows, audits = [], [], [], [], []
    for config in configs:
        targets = {}
        for acc in DATASETS:
            paths = [Path(str(stems(acc, s, config, count_source)) + '__target.tsv.gz') for s in samples(acc)]
            if not all(p.exists() for p in paths):
                return False
            targets[acc] = pd.concat([pd.read_csv(p, sep='\t') for p in paths], ignore_index=True)
            for s in samples(acc):
                ap = Path(str(stems(acc, s, config, count_source)) + '.json')
                audits.append(json.loads(ap.read_text(encoding='utf-8')))
        common = None
        scores, manifests = {}, {}
        for acc in DATASETS:
            m = targets[acc][['sample', 'condition']].drop_duplicates().sort_values('sample')
            manifests[acc] = m
            scores[acc] = {}
            for metric in METRICS:
                mat = targets[acc].pivot(index=KEY, columns='sample', values=metric).reindex(columns=m['sample'])
                complete = mat.index[np.isfinite(mat.to_numpy()).all(axis=1)]
                common = complete if common is None else common.intersection(complete)
                scores[acc][metric] = mat
        common = common.sort_values()
        for metric in METRICS:
            effs = []
            for acc in DATASETS:
                mat = scores[acc][metric].loc[common]
                case = manifests[acc].condition.to_numpy() == 'MCAO'
                x = mat.to_numpy()
                eff = x[:, case].mean(axis=1) - x[:, ~case].mean(axis=1)
                eff[np.abs(eff) <= TOL] = 0; effs.append(eff)
                row = common.to_frame(index=False)
                row.insert(0, 'metric', metric); row.insert(0, 'dataset', acc)
                row.insert(0, 'config', config); row.insert(0, 'count_source', count_source)
                row['mcao_minus_sham'] = eff; effects.append(row.rename(columns=RENAMES))
                long = mat.rename_axis(columns='sample').stack().rename('value').reset_index().rename(columns=RENAMES)
                long.insert(0, 'metric', metric); long.insert(0, 'dataset', acc)
                long.insert(0, 'config', config); long.insert(0, 'count_source', count_source)
                means.append(long.merge(manifests[acc], on='sample', validate='many_to_one'))
            a, b = effs
            nz = (a != 0) & (b != 0); same = nz & (np.sign(a) == np.sign(b))
            all_rows.append({'count_source': count_source, 'config': config, 'metric': metric,
                'n_candidates': len(common), 'n_nonzero_both': int(nz.sum()), 'n_same_direction': int(same.sum()),
                'n_zero_discovery': int((a == 0).sum()), 'n_zero_validation': int((b == 0).sum()),
                'same_direction_all_fraction': float(same.mean()),
                'same_direction_nonzero_fraction': float(same.sum() / nz.sum()),
                'spearman_rho': float(stats.spearmanr(a, b).statistic)})
        for acc in DATASETS:
            d = targets[acc]
            for sample, b in d.groupby('sample', sort=False):
                for universe, t in [('all_reported_target_edges', b), ('fixed_common_target_edges', b.set_index(KEY).loc[common].reset_index())]:
                    for metric in ['magnitude_rank', 'diagnostic_unique_column_magnitude_rank']:
                        r = t[metric].to_numpy()
                        counts = np.unique(r, return_counts=True)[1]
                        tie_rows.append({'count_source': count_source, 'config': config, 'dataset': acc,
                            'sample': sample, 'universe': universe, 'metric': metric, 'n_edges': len(r),
                            'n_rank_one': int((np.abs(r - 1) <= TOL).sum()),
                            'rank_one_fraction': float((np.abs(r - 1) <= TOL).mean()),
                            'n_edges_in_exact_tie_groups': int(counts[counts > 1].sum()),
                            'largest_exact_tie_group': int(counts.max())})
    TABLES.mkdir(parents=True, exist_ok=True)
    prefix = 'rank_diagnostic_' + count_source + '_'
    for suffix, frame in [('summary.tsv', pd.DataFrame(all_rows)), ('effects.tsv.gz', pd.concat(effects, ignore_index=True)),
                          ('animal_scores.tsv.gz', pd.concat(means, ignore_index=True)), ('ties.tsv', pd.DataFrame(tie_rows))]:
        frame.to_csv(TABLES / (prefix + suffix), sep='\t', index=False)
    meta = {'status': 'complete_for_requested_configs', 'count_source': count_source, 'configs': configs,
            'n_library_config_audits': len(audits), 'native_raw_reruns_match_prior': all(
                a['raw_target_scores_verified_against_prior'] for a in audits) if count_source == 'raw' else None,
            'full_network_intervention': 'each unique magnitude score column ranked once, on the identical complete native network',
            'inferential_pvalues_computed': False, 'audits': audits}
    (TABLES / (prefix + 'audit.json')).write_text(json.dumps(meta, indent=2, allow_nan=False), encoding='utf-8')
    print(pd.DataFrame(all_rows).to_string(index=False), flush=True)
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count-source', choices=['raw', 'decontx'], default='raw')
    parser.add_argument('--datasets', nargs='+', choices=DATASETS, default=DATASETS)
    parser.add_argument('--configs', nargs='+', choices=['primary', 'reference_singlet'], default=['primary', 'reference_singlet'])
    parser.add_argument('--samples', nargs='+')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    completed = []
    for config in args.configs:
        for acc in args.datasets:
            for sample in (args.samples or samples(acc)):
                process(acc, sample, config, args.count_source, args.force)
        completed.append(config)
        summarize(args.count_source, completed)
