"""Isolated fixed-commit upstream aggregation audit on saved small network tables."""
from __future__ import annotations
import ast
import gc
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import sys
for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[name] = '1'
import numpy as np
import pandas as pd
import psutil
from scipy import stats

WORK = Path(__file__).resolve().parent
LIT = WORK / 'repro_literature'
TABLES = WORK.parent / 'outputs/reproducibility_v2/tables'
KEY = ['source', 'target', 'ligand_complex', 'receptor_complex']
RENAMES = dict(zip(KEY, ['sender', 'receiver', 'ligand', 'receptor']))
SCORES = ['lr_means', 'expr_prod', 'lrscore']
FUNCS = ['_rank_aggregate', '_corr_beta_pvals', '_rho_scores', '_robust_rank_aggregate']
TOL = 1e-12
CAP = 400 * 1024 ** 2


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def memory():
    x = psutil.Process().memory_info()
    peak = getattr(x, 'peak_wset', x.rss)
    if peak > CAP:
        raise RuntimeError('Isolated audit exceeded its 400 MiB memory bound')
    return {'rss_bytes': x.rss, 'peak_bytes': peak}


def load_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    provenance_path = LIT / 'liana_pr261_source_provenance.json'
    provenance = json.loads(provenance_path.read_text(encoding='utf-8'))
    assert provenance['merged'] and provenance['merge_commit_sha'] == 'd4211373692e7b9c10210488ccb1efe06452b097'
    for rec in provenance['sources']:
        assert sha(LIT / rec['filename']) == rec['sha256']
    source = LIT / 'liana_merge_d421137_aggregate.py'
    tree = ast.parse(source.read_text(encoding='utf-8'))
    funcs = [x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name in FUNCS]
    assert [x.name for x in funcs] == FUNCS
    header = '''# Derived executable subset of scverse/liana merge d4211373692e7b9c10210488ccb1efe06452b097.
# Function bodies are AST-identical to the archived upstream source. No new LIANA package is loaded.
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import beta, rankdata
def _logg(*args, **kwargs):
    raise RuntimeError('Logging branch was not expected for this three-column aggregation audit')
'''
    extracted = LIT / 'liana_merge_d421137_extracted_aggregation.py'
    extracted.write_text(header + '\n\n'.join(ast.unparse(x) for x in funcs) + '\n', encoding='utf-8')
    derived = ast.parse(extracted.read_text(encoding='utf-8'))
    extracted_funcs = {x.name: x for x in derived.body if isinstance(x, ast.FunctionDef)}
    assert all(ast.dump(x, include_attributes=False) == ast.dump(extracted_funcs[x.name], include_attributes=False) for x in funcs)
    upstream = load_file('fixed_commit_upstream_subset', extracted)
    old_path = Path(importlib.util.find_spec('liana').origin).parent / 'method/_pipe_utils/_aggregate.py'
    old = load_file('installed_110_aggregation_only', old_path)
    keys = pd.read_csv(TABLES / 'null_fixed_candidate_keys.tsv', sep='\t').rename(columns={v: k for k, v in RENAMES.items()})
    inputs = {str(source.relative_to(WORK)): sha(source), str(extracted.relative_to(WORK)): sha(extracted),
              str(provenance_path.relative_to(WORK)): sha(provenance_path),
              'repro_upstream_correction_plan.md': sha(WORK / 'repro_upstream_correction_plan.md'),
              str(old_path): sha(old_path)}
    rows, animal, memory_checks = [], [], [memory()]
    for config, n in [('primary', 187), ('reference_singlet', 174)]:
        target = pd.MultiIndex.from_frame(keys[keys.config == config][KEY]).sort_values()
        assert len(target) == n
        for acc, expected in [('GSE174574', 6), ('GSE245386', 5)]:
            paths = sorted((WORK / 'repro_liana_rank_diagnostic/raw' / acc).glob(config + '__*__full_network.tsv.gz'))
            assert len(paths) == expected
            for path in paths:
                d = pd.read_csv(path, sep='\t', float_precision='round_trip',
                                usecols=KEY + SCORES + ['sample', 'condition', 'magnitude_rank',
                                                       'diagnostic_unique_column_magnitude_rank']).set_index(KEY)
                assert not d.index.duplicated().any()
                d[SCORES] = d[SCORES].astype(np.float32)
                assert np.isfinite(d[SCORES]).all().all()
                sample = d['sample'].iloc[0]
                ap = path.parent / (config + '__' + sample + '.json')
                a = json.loads(ap.read_text(encoding='utf-8'))
                assert sha(path) == a['full_network_sha256']
                inputs[str(path.relative_to(WORK))] = sha(path)
                inputs[str(ap.relative_to(WORK))] = sha(ap)
                native_specs, unique_specs = a['native_magnitude_specs'], a['diagnostic_unique_score_specs']
                legacy_native = old._rank_aggregate(d.reset_index().copy(), native_specs, 'rra')
                legacy_unique = old._rank_aggregate(d.reset_index().copy(), unique_specs, 'rra')
                assert np.max(np.abs(legacy_native - d.magnitude_rank)) <= TOL
                assert np.max(np.abs(legacy_unique - d.diagnostic_unique_column_magnitude_rank)) <= TOL
                before = d[SCORES].copy()
                fixed = upstream._rank_aggregate(d.reset_index(), native_specs, 'rra')
                pd.testing.assert_frame_equal(before, d[SCORES])
                # Same restored float32 score values; float64 only changes the
                # destination dtype for the old function's in-place ranks.
                wide = d.reset_index().copy()
                wide[SCORES] = wide[SCORES].astype(np.float64)
                legacy_unique_float64 = old._rank_aggregate(wide.copy(), unique_specs, 'rra')
                upstream_float64 = upstream._rank_aggregate(wide.copy(), native_specs, 'rra')
                rank_exact = all(np.array_equal(stats.rankdata(d[col], method='average'),
                                                stats.rankdata(wide[col], method='average')) for col in SCORES)
                assert rank_exact
                assert np.array_equal(legacy_unique_float64, upstream_float64)
                positions = d.index.get_indexer(target)
                assert (positions >= 0).all()
                for scope, idx in [('complete_full_network', np.arange(len(d))), ('fixed_target_set', positions)]:
                    diff = np.abs(fixed[idx] - legacy_unique[idx])
                    rows.append({'config': config, 'dataset': acc, 'sample': sample, 'scope': scope,
                                 'n_edges': len(idx), 'unique_rank_vectors_equal': rank_exact,
                                 'upstream_vs_legacy_unique_bitwise_equal': bool(np.array_equal(fixed[idx], legacy_unique[idx])),
                                 'max_abs_upstream_vs_legacy_unique': float(diff.max()),
                                 'n_difference_gt_1e12': int((diff > 1e-12).sum()),
                                 'n_difference_gt_1e6': int((diff > 1e-6).sum()),
                                 'float32_rankdata_dtype': str(stats.rankdata(d[SCORES[0]], method='average').dtype),
                                 'float64_rankdata_dtype': str(stats.rankdata(wide[SCORES[0]], method='average').dtype),
                                 'legacy_float64_vs_upstream_float64_bitwise_equal': bool(np.array_equal(upstream_float64[idx], legacy_unique_float64[idx])),
                                 'max_abs_upstream_float32_vs_float64': float(np.abs(fixed[idx] - upstream_float64[idx]).max())})
                for metric, values in [('native_110_priority', 1 - legacy_native[positions]),
                                       ('legacy_unique_column_priority', 1 - legacy_unique[positions]),
                                       ('fixed_merge_upstream_priority', 1 - fixed[positions])]:
                    r = target.to_frame(index=False).rename(columns=RENAMES)
                    r['config'], r['dataset'], r['sample'], r['condition'], r['metric'] = config, acc, sample, a['condition'], metric
                    r['value'] = values
                    animal.append(r)
                memory_checks.append(memory())
                print('AUDITED', config, sample, 'edges', len(d), 'max delta', float(np.abs(fixed - legacy_unique).max()), flush=True)
                del d, wide, before, legacy_native, legacy_unique, fixed, legacy_unique_float64, upstream_float64
                gc.collect()
    animals = pd.concat(animal, ignore_index=True)
    outkey = list(RENAMES.values())
    effects, summaries = [], []
    for (config, metric), block in animals.groupby(['config', 'metric'], sort=False):
        deltas = []
        for acc in ['GSE174574', 'GSE245386']:
            d = block[block.dataset == acc]
            m = d[['sample', 'condition']].drop_duplicates().sort_values('sample')
            assert m.condition.value_counts().to_dict() == {'Sham': 3, 'MCAO': 3 if acc == 'GSE174574' else 2}
            x = d.pivot(index=outkey, columns='sample', values='value').reindex(columns=m['sample'])
            case = m.condition.to_numpy() == 'MCAO'
            delta = x.to_numpy()[:, case].mean(axis=1) - x.to_numpy()[:, ~case].mean(axis=1)
            delta[np.abs(delta) <= TOL] = 0
            deltas.append(delta)
            row = x.index.to_frame(index=False)
            row['config'], row['dataset'], row['metric'], row['mcao_minus_sham'] = config, acc, metric, delta
            effects.append(row)
        av, bv = deltas
        nz = (av != 0) & (bv != 0)
        same = nz & (np.sign(av) == np.sign(bv))
        summaries.append({'config': config, 'metric': metric, 'n_candidates': len(av),
                          'n_nonzero_both': int(nz.sum()), 'n_same_direction': int(same.sum()),
                          'same_direction_all_fraction': float(same.mean()),
                          'same_direction_nonzero_fraction': float(same.sum() / nz.sum()),
                          'n_zero_discovery': int((av == 0).sum()), 'n_zero_validation': int((bv == 0).sum()),
                          'spearman_rho': float(stats.spearmanr(av, bv).statistic)})
    comparisons, summary = pd.DataFrame(rows), pd.DataFrame(summaries)
    tables = {'per_library_equality.tsv': comparisons, 'animal_scores.tsv.gz': animals,
              'effects.tsv.gz': pd.concat(effects, ignore_index=True), 'summary.tsv': summary}
    for suffix, table in tables.items():
        table.to_csv(TABLES / ('rank_upstream_' + suffix), sep='\t', index=False)
    audit = {'status': 'complete', 'n_library_config_audits': 22, 'upstream_pr': provenance['pr_url'],
             'upstream_merged_at': provenance['merged_at'], 'upstream_merge_commit': provenance['merge_commit_sha'],
             'function_bodies_AST_identical': True, 'extracted_functions': FUNCS,
             'new_package_installed': False, 'installed_package_modified': False, 'new_full_pipeline_run': False,
             'aggregation_only_on_unchanged_saved_networks': True, 'score_dtype_restored': 'float32',
             'known_upstream_fix_not_a_new_bug_discovery': True, 'new_pvalues': False,
             'all_unique_rank_vectors_equal': bool(comparisons.unique_rank_vectors_equal.all()),
             'all_legacy_float64_vs_upstream_float64_bitwise_equal': bool(comparisons.legacy_float64_vs_upstream_float64_bitwise_equal.all()),
             'all_upstream_vs_original_unique_diagnostic_bitwise_equal': bool(comparisons.upstream_vs_legacy_unique_bitwise_equal.all()),
             'rankdata_float32_runtime_dtype': 'float32',
             'maximum_absolute_float32_vs_float64_sensitivity': float(comparisons.max_abs_upstream_float32_vs_float64.max()),
             'maximum_absolute_difference_from_original_unique_diagnostic': float(comparisons.max_abs_upstream_vs_legacy_unique.max()),
             'all_differences_below_1e6': bool(comparisons.n_difference_gt_1e6.eq(0).all()),
             'memory_cap_bytes': CAP, 'maximum_observed_process_peak_bytes': max(x['peak_bytes'] for x in memory_checks),
             'source_sha256': inputs, 'script_sha256': sha(Path(__file__)),
             'output_sha256': {suffix: sha(TABLES / ('rank_upstream_' + suffix)) for suffix in tables},
             'versions': {x: importlib.metadata.version(x) for x in ['numpy', 'pandas', 'scipy', 'liana']},
             'limits': ['The merge also changes defaults and other pipeline code; those changes were not executed.',
                        'Numerical equality after harmonizing rank-matrix precision is not a claim of full-version pipeline equivalence.']}
    (TABLES / 'rank_upstream_audit.json').write_text(json.dumps(audit, indent=2, allow_nan=False), encoding='utf-8')
    print(summary.to_string(index=False), flush=True)
    print(json.dumps({k: audit[k] for k in ['maximum_absolute_difference_from_original_unique_diagnostic',
                     'all_upstream_vs_original_unique_diagnostic_bitwise_equal',
                     'all_legacy_float64_vs_upstream_float64_bitwise_equal', 'maximum_observed_process_peak_bytes']}, indent=2))


if __name__ == '__main__':
    main()
