"""Descriptive fixed-universe rerank/context audit; no additional P values."""
from __future__ import annotations
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
for _n in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_n] = '1'
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'outputs/reproducibility_v2/tables'
KEY = ['sender', 'receiver', 'ligand', 'receptor']
DATASETS = ['GSE174574', 'GSE245386']
CONFIGS = ['primary', 'reference_singlet']
EXPR = ['custom_score', 'lr_means', 'expr_prod', 'lrscore']
TOL = 1e-12
INPUTS = {}
CHECKS = []

def register(p):
    name = str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)
    INPUTS[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p

def read(p, **kwargs):
    return pd.read_csv(register(p), sep='\t', **kwargs)

def check(name, condition):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)

def contrast(x, conditions):
    case = np.asarray(conditions) == 'MCAO'
    a = x[:, case].mean(axis=1) - x[:, ~case].mean(axis=1)
    a[np.abs(a) <= TOL] = 0
    return a

def concordance(a, b):
    nz = (a != 0) & (b != 0)
    same = nz & (np.sign(a) == np.sign(b))
    return {'n_candidates': len(a), 'n_nonzero_both': int(nz.sum()),
            'n_same_direction': int(same.sum()), 'n_zero_discovery': int((a == 0).sum()),
            'n_zero_validation': int((b == 0).sum()),
            'same_direction_all_fraction': float(same.sum() / len(a)),
            'same_direction_nonzero_fraction': float(same.sum() / nz.sum()),
            'spearman_rho': float(stats.spearmanr(a, b).statistic)}

def tie_summary(rank):
    rank = np.asarray(rank, float)
    _, counts = np.unique(rank, return_counts=True)
    counts_rounded = np.unique(np.round(rank, 12), return_counts=True)[1]
    return {'n_target_edges': len(rank), 'n_unique_rank_values': len(counts),
            'n_edges_in_exact_tie_groups': int(counts[counts > 1].sum()),
            'fraction_edges_in_exact_tie_groups': float(counts[counts > 1].sum() / len(rank)),
            'largest_exact_tie_group': int(counts.max()),
            'n_edges_in_12decimal_tie_groups': int(counts_rounded[counts_rounded > 1].sum()),
            'n_native_rank_one': int((np.abs(rank - 1) <= TOL).sum()),
            'native_rank_one_fraction': float((np.abs(rank - 1) <= TOL).mean()),
            'n_native_rank_zero': int((np.abs(rank) <= TOL).sum()),
            'native_rank_zero_fraction': float((np.abs(rank) <= TOL).mean()),
            'native_rank_min': float(rank.min()), 'native_rank_max': float(rank.max())}

def implementation_probe():
    pkg = Path(importlib.util.find_spec('liana').origin).parent
    init = ast.parse(register(pkg / 'method/__init__.py').read_text(encoding='utf-8'))
    default_methods = next([n.id for n in node.value.elts] for node in init.body
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == '_methods' for t in node.targets))
    specs = {}
    for method in default_methods:
        tree = ast.parse(register(pkg / 'method/sc' / ('_' + method + '.py')).read_text(encoding='utf-8'))
        call = next(n for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                    and n.func.id == 'MethodMeta')
        kw = {k.arg: k.value for k in call.keywords}
        magnitude = ast.literal_eval(kw['magnitude'])
        if magnitude is not None:
            specs[ast.literal_eval(kw['method_name'])] = (magnitude, ast.literal_eval(kw['magnitude_ascending']))
    check('default_shared_expr_prod_occurs_twice', sum(v[0] == 'expr_prod' for v in specs.values()) == 2)
    unique_specs = {}
    used = set()
    for k, value in specs.items():
        if value[0] not in used:
            unique_specs[k] = value; used.add(value[0])
    path = register(pkg / 'method/_pipe_utils/_aggregate.py')
    loader = importlib.util.spec_from_file_location('official_liana_rank_audit', path)
    module = importlib.util.module_from_spec(loader); loader.loader.exec_module(module)
    toy = pd.DataFrame({'lr_means': [1., 2., 3., 4.], 'expr_prod': [1., 4., 9., 16.],
                        'lrscore': [.1, .2, .3, .4]})
    native = toy.copy(); unique = toy.copy()
    native_rra = module._rank_aggregate(native, specs, 'rra')
    unique_rra = module._rank_aggregate(unique, unique_specs, 'rra')
    check('native_probe_shared_column_reverses_descending_ranks', np.array_equal(native.expr_prod, [1, 2, 3, 4]))
    check('diagnostic_unique_column_has_expected_descending_ranks', np.array_equal(unique.expr_prod, [4, 3, 2, 1]))
    result = toy.add_prefix('input_')
    result = pd.concat([result, native.add_prefix('native_intermediate_rank_'),
                        unique.add_prefix('unique_column_diagnostic_rank_')], axis=1)
    result['native_rra'] = native_rra
    result['unique_column_diagnostic_rra'] = unique_rra
    result.insert(0, 'toy_row', np.arange(1, len(result) + 1))
    info = {'default_methods': default_methods, 'native_magnitude_specs': specs,
            'unique_column_diagnostic_specs': unique_specs,
            'official_tagged_source': 'https://raw.githubusercontent.com/saezlab/liana-py/v1.10.0/src/liana/method/_pipe_utils/_aggregate.py',
            'native_source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'package_modified': False, 'existing_native_scores_replaced': False,
            'diagnostic_is_real_network_rerun': False,
            'interpretation': 'The shared expr_prod column is transformed twice in the official default magnitude loop, reversing its first descending rank. Real-network quantitative impact is not established by this four-row probe.'}
    return result, info

def main():
    register(ROOT / 'repro_rank_mechanism_plan.md')
    src = read(OUT / 'null_fixed_common_animal_scores.tsv.gz')
    keys = read(OUT / 'null_fixed_candidate_keys.tsv')
    summaries, effects, animal_values = [], [], []
    all_context, all_types, all_ties, all_zeros = [], [], [], []
    for config in CONFIGS:
        fixed = pd.MultiIndex.from_frame(keys[keys.config == config][KEY]).sort_values()
        by_metric = {}
        for metric in EXPR + ['magnitude_priority']:
            for representation in (['native_expression', 'fixed_universe_percentile'] if metric in EXPR else ['native_liana_consensus']):
                cohort_effects = []
                for acc in DATASETS:
                    d = src[(src.config == config) & (src.dataset == acc) & (src.metric == metric)]
                    m = d[['sample', 'condition']].drop_duplicates().sort_values('sample')
                    mat = d.pivot(index=KEY, columns='sample', values='value').reindex(index=fixed, columns=m['sample'])
                    check(f'{config}/{acc}/{metric}/{representation}/finite_complete', np.isfinite(mat).all().all())
                    if representation == 'fixed_universe_percentile':
                        mat = mat.rank(axis=0, method='average', pct=True)
                        saved = d.pivot(index=KEY, columns='sample', values='within_fixed_universe_percentile').reindex(index=fixed, columns=m['sample'])
                        check(f'{config}/{acc}/{metric}/rerank_matches_saved_interface', np.allclose(mat, saved, atol=TOL, rtol=0))
                    eff = contrast(mat.to_numpy(float), m.condition.to_numpy())
                    cohort_effects.append(eff)
                    row = fixed.to_frame(index=False)
                    row.insert(0, 'representation', representation); row.insert(0, 'metric', metric)
                    row.insert(0, 'dataset', acc); row.insert(0, 'config', config)
                    row['mcao_minus_sham'] = eff; effects.append(row)
                    values = mat.rename_axis(columns='sample').stack().rename('value').reset_index()
                    values.insert(0, 'representation', representation); values.insert(0, 'metric', metric)
                    values.insert(0, 'dataset', acc); values.insert(0, 'config', config)
                    animal_values.append(values.merge(m, on='sample', validate='many_to_one'))
                summary = {'config': config, 'metric': metric, 'representation': representation} | concordance(*cohort_effects)
                summaries.append(summary)
                by_metric[metric, representation] = summary
        for acc in DATASETS:
            d = src[(src.config == config) & (src.dataset == acc) & (src.metric == 'magnitude_priority')]
            m = d[['sample', 'condition']].drop_duplicates().sort_values('sample')
            fixed_native_rank = []
            for sample, condition in m.itertuples(index=False, name=None):
                stem = ROOT / 'revision_results/liana' / acc / (config + '__' + sample)
                ap = register(Path(str(stem) + '.json'))
                audit = json.loads(ap.read_text(encoding='utf-8'))
                t = read(Path(str(stem) + '.tsv.gz')).rename(columns={'source': 'sender', 'target': 'receiver',
                    'ligand_complex': 'ligand', 'receptor_complex': 'receptor'}).set_index(KEY)
                meta = {'config': config, 'dataset': acc, 'sample': sample, 'condition': condition}
                check(f'{config}/{acc}/{sample}/target_export_size_matches_audit', len(t) == audit['exported_target_rows'])
                all_context.append(meta | {'n_eligible_full_network_rows': audit['all_pair_rows'],
                    'n_reported_target_direction_rows': len(t), 'n_context_cells': audit['n_context_cells'],
                    'n_context_types': len(audit['context_counts']),
                    'context_type_names': '|'.join(sorted(audit['context_counts']))})
                for typ, count in audit['context_counts'].items():
                    all_types.append(meta | {'context_cell_type': typ, 'n_cells': count})
                for universe, rank in [('all_reported_target_edges', t.magnitude_rank.to_numpy()),
                                       ('fixed_common_target_edges', t.loc[fixed].magnitude_rank.to_numpy())]:
                    all_ties.append(meta | {'universe': universe} | tie_summary(rank))
                fixed_native_rank.append(t.loc[fixed].magnitude_rank.to_numpy())
            rank = np.column_stack(fixed_native_rank)
            native_priority = 1 - rank
            effect = contrast(native_priority, m.condition.to_numpy())
            all_at_one = (np.abs(rank - 1) <= TOL).all(axis=1)
            any_at_one = (np.abs(rank - 1) <= TOL).any(axis=1)
            all_zeros.append({'config': config, 'dataset': acc, 'n_fixed_candidates': len(fixed),
                'n_zero_priority_disease_effect': int((effect == 0).sum()),
                'n_edges_at_rank_one_in_all_animals': int(all_at_one.sum()),
                'n_edges_at_rank_one_in_any_animal': int(any_at_one.sum()),
                'fraction_edges_at_rank_one_in_all_animals': float(all_at_one.mean()),
                'n_zero_effect_explained_by_all_animal_rank_one': int(((effect == 0) & all_at_one).sum())})
    comparison = pd.DataFrame(summaries)
    context = pd.DataFrame(all_context)
    ties = pd.DataFrame(all_ties)
    ranges = []
    for (config, acc), b in context.groupby(['config', 'dataset'], sort=False):
        rec = {'config': config, 'dataset': acc, 'n_animals': len(b),
               'n_distinct_context_type_sets': int(b.context_type_names.nunique())}
        for col in ['n_eligible_full_network_rows', 'n_reported_target_direction_rows', 'n_context_cells', 'n_context_types']:
            rec[col + '_min'] = int(b[col].min()); rec[col + '_max'] = int(b[col].max())
        for universe in ['all_reported_target_edges', 'fixed_common_target_edges']:
            r = ties[(ties.config == config) & (ties.dataset == acc) & (ties.universe == universe)]
            for col in ['native_rank_one_fraction', 'fraction_edges_in_exact_tie_groups']:
                rec[universe + '__' + col + '_min'] = float(r[col].min())
                rec[universe + '__' + col + '_max'] = float(r[col].max())
        ranges.append(rec)
    probe, probe_info = implementation_probe()
    tables = {'comparison_summary.tsv': comparison, 'figure_source_comparison.tsv': comparison.copy(),
              'observed_effects.tsv.gz': pd.concat(effects, ignore_index=True),
              'per_animal_representations.tsv.gz': pd.concat(animal_values, ignore_index=True),
              'network_context_per_animal.tsv': context, 'network_context_cell_counts.tsv': pd.DataFrame(all_types),
              'native_target_ties_per_animal.tsv': ties, 'native_zero_effects.tsv': pd.DataFrame(all_zeros),
              'network_context_ranges.tsv': pd.DataFrame(ranges),
              'implementation_four_row_probe.tsv': probe}
    for name, table in tables.items():
        table.to_csv(OUT / ('rank_' + name), sep='\t', index=False)
    meta = {'status': 'complete', 'n_checks': len(CHECKS), 'checks': CHECKS, 'input_sha256': INPUTS,
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'pvalues_computed': False, 'fixed_percentiles_recompute_full_network_consensus': False,
            'target_tie_distributions_represent_full_network': False,
            'implementation_probe': probe_info, 'table_rows': {n: len(d) for n, d in tables.items()}}
    (OUT / 'rank_audit.json').write_text(json.dumps(meta, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    (OUT / 'rank_methods.md').write_text((ROOT / 'repro_rank_mechanism_plan.md').read_text(encoding='utf-8'), encoding='utf-8')
    print(comparison.to_string(index=False)); print(pd.DataFrame(ranges).to_string(index=False)); print(pd.DataFrame(all_zeros).to_string(index=False))

if __name__ == '__main__':
    main()
