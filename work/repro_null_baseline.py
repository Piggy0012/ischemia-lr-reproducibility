"""Exhaustive animal-label nulls preserving whole LR vectors and shared genes.

Run with work/revision_env/Scripts/python.exe work/repro_null_baseline.py.
The amendment is in repro_null_plan.md and must precede these outputs.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path

for _name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_name] = '1'

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'outputs/reproducibility_v2/tables'
DATASETS = ['GSE174574', 'GSE245386']
CONFIGS = ['primary', 'reference_singlet']
KEY = ['sender', 'receiver', 'ligand', 'receptor']
LMETRICS = ['lr_means', 'expr_prod', 'lrscore', 'magnitude_priority']
METRICS = ['custom_score'] + LMETRICS
EXPR = METRICS[:-1]
STATISTICS = ['same_direction_all_fraction', 'same_direction_nonzero_fraction', 'spearman_rho']
TOL = 1e-12
INPUTS = {}
CHECKS = []


def register(path):
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read(path, **kwargs):
    return pd.read_csv(register(path), sep='\t', **kwargs)


def check(name, value, detail=None):
    CHECKS.append({'name': name, 'passed': bool(value), 'detail': detail})
    if not value:
        raise AssertionError(name + (': ' + str(detail) if detail else ''))


def manifest(acc, config):
    folder = ROOT / ('processed' if config == 'primary' else 'revision_analysis/processed') / acc
    d = read(folder / 'sample_cell_counts.tsv')
    d = d.loc[d.config == config, ['sample', 'condition']].drop_duplicates().sort_values('sample')
    check(f'{acc}/{config}/unique_animals', not d['sample'].duplicated().any())
    check(f'{acc}/{config}/expected_group_counts',
          d.condition.value_counts().to_dict() == {'Sham': 3, 'MCAO': 3 if acc == DATASETS[0] else 2})
    return folder, d.reset_index(drop=True)


def all_custom(acc, config, folder, m, lr):
    """Build all resource- and gene-covered edges before any disease group gate."""
    sample_order = m['sample'].tolist()
    ms, fs = [], []
    usecols = [config + '__' + typ for typ in ['Astrocyte', 'Endothelial']]
    for sample in sample_order:
        a = read(folder / (sample + '_mean_logcp10k.tsv.gz'), index_col=0)[usecols]
        b = read(folder / (sample + '_fractions.tsv.gz'), index_col=0)[usecols]
        check(f'{acc}/{config}/{sample}/matching_mean_fraction_genes', a.index.equals(b.index))
        ms.append(a); fs.append(b)
    genes = ms[0].index
    check(f'{acc}/{config}/matching_animal_gene_order', all(x.index.equals(genes) for x in ms + fs))
    pos = {g: j for j, g in enumerate(genes)}
    means = {typ: np.column_stack([x[config + '__' + typ].to_numpy(float) for x in ms])
             for typ in ['Astrocyte', 'Endothelial']}
    fractions = {typ: np.column_stack([x[config + '__' + typ].to_numpy(float) for x in fs])
                 for typ in ['Astrocyte', 'Endothelial']}
    keys, scores, ldet, rdet = [], [], [], []
    for sender, receiver in [('Astrocyte', 'Endothelial'), ('Endothelial', 'Astrocyte')]:
        for ligand, receptor in lr[['source_genesymbol', 'target_genesymbol']].itertuples(index=False, name=None):
            lg, rg = ligand.split('_'), receptor.split('_')
            if not all(g in pos for g in lg + rg):
                continue
            li, ri = [pos[g] for g in lg], [pos[g] for g in rg]
            keys.append((sender, receiver, ligand, receptor))
            scores.append(np.sqrt(means[sender][li].min(axis=0) * means[receiver][ri].min(axis=0)))
            ldet.append(fractions[sender][li].min(axis=0))
            rdet.append(fractions[receiver][ri].min(axis=0))
    index = pd.MultiIndex.from_tuples(keys, names=KEY)
    result = {'scores': pd.DataFrame(np.array(scores), index=index, columns=sample_order),
              'ligand_detection': np.array(ldet), 'receptor_detection': np.array(rdet),
              'manifest': m}
    check(f'{acc}/{config}/finite_all_resource_scores', np.isfinite(result['scores']).all().all())
    # Independently compare reconstructed scores to the previously exported table.
    base = ROOT / ('results' if config == 'primary' else 'revision_analysis/results') / acc
    prior = read(base / (config + '__communication_sample_scores.tsv.gz'))
    prior_matrix = prior.pivot(index=KEY, columns='sample', values='score').reindex(columns=sample_order)
    error = np.max(np.abs(result['scores'].loc[prior_matrix.index].to_numpy() - prior_matrix.to_numpy()))
    check(f'{acc}/{config}/reconstructed_custom_matches_prior', error < TOL, float(error))
    return result


def liana(acc, config, m):
    rows = []
    for sample, condition in m.itertuples(index=False, name=None):
        stem = ROOT / 'revision_results/liana' / acc / (config + '__' + sample)
        audit = json.loads(register(Path(str(stem) + '.json')).read_text(encoding='utf-8'))
        check(f'{acc}/{config}/{sample}/liana_full_context_eligible_only',
              audit['computed_all_cell_type_pairs'] is True and audit['parameters']['return_all_lrs'] is False)
        d = read(Path(str(stem) + '.tsv.gz')).rename(columns={
            'source': 'sender', 'target': 'receiver', 'ligand_complex': 'ligand', 'receptor_complex': 'receptor'})
        check(f'{acc}/{config}/{sample}/liana_unique_correct_sample',
              not d.duplicated(KEY).any() and set(d['sample']) == {sample} and set(d.condition) == {condition})
        d['magnitude_priority'] = 1 - d.magnitude_rank
        rows.append(d)
    all_rows = pd.concat(rows, ignore_index=True)
    matrices = {metric: all_rows.pivot(index=KEY, columns='sample', values=metric).reindex(columns=m['sample'])
                for metric in LMETRICS}
    return matrices


def allocations(m):
    ncase = int((m.condition == 'MCAO').sum())
    result = []
    observed = (m.condition == 'MCAO').to_numpy()
    for idx, selected in enumerate(itertools.combinations(range(len(m)), ncase)):
        mask = np.zeros(len(m), bool); mask[list(selected)] = True
        result.append({'index': idx, 'mask': mask, 'is_observed': np.array_equal(mask, observed),
                       'mcao_animals': '|'.join(m.loc[mask, 'sample'])})
    return result


def delta(matrix, case):
    x = matrix.to_numpy(float) if hasattr(matrix, 'to_numpy') else matrix
    result = x[:, case].mean(axis=1) - x[:, ~case].mean(axis=1)
    result[np.abs(result) <= TOL] = 0
    return result


def comparison(a, b):
    check_finite = np.isfinite(a).all() and np.isfinite(b).all()
    if not check_finite or not len(a):
        raise ValueError('No missing or empty effect vectors are allowed')
    nonzero = (a != 0) & (b != 0)
    same = nonzero & (np.sign(a) == np.sign(b))
    rho = float(stats.spearmanr(a, b).statistic)
    if not nonzero.any() or not np.isfinite(rho):
        raise ValueError('Undefined planned statistic')
    return {'n_candidates': len(a), 'n_nonzero_both': int(nonzero.sum()),
            'n_zero_discovery': int((a == 0).sum()), 'n_zero_validation': int((b == 0).sum()),
            'n_same_direction': int(same.sum()),
            'same_direction_all_fraction': float(same.sum() / len(a)),
            'same_direction_nonzero_fraction': float(same.sum() / nonzero.sum()),
            'spearman_rho': rho}


def group_gate(data, mask, cutoff=.1):
    ld, rd = data['ligand_detection'] >= cutoff, data['receptor_detection'] >= cutoff
    lp = np.maximum(ld[:, mask].sum(axis=1), ld[:, ~mask].sum(axis=1)) >= 2
    rp = np.maximum(rd[:, mask].sum(axis=1), rd[:, ~mask].sum(axis=1)) >= 2
    return lp & rp


def summarize(long, groupcols, family):
    result = []
    for key, block in long.groupby(groupcols, sort=False, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        obs = block.loc[block.is_observed]
        if len(obs) != 1 or len(block) != 200:
            raise AssertionError('Need exactly 200 allocations and exactly one observed allocation')
        vals = block['value'].to_numpy(float)
        if not np.isfinite(vals).all():
            raise AssertionError('An intended null statistic is undefined')
        observed = float(obs['value'].iloc[0])
        p = float((vals >= observed - TOL).mean())
        result.append(dict(zip(groupcols, key)) | {
            'test_family': family, 'observed': observed, 'null_mean': float(vals.mean()),
            'null_median': float(np.median(vals)), 'null_q025': float(np.quantile(vals, .025)),
            'null_q975': float(np.quantile(vals, .975)), 'n_allocations': len(vals),
            'n_greater_or_equal_observed': int((vals >= observed - TOL).sum()),
            'exact_upper_tail_p': p, 'observed_n_candidates': int(obs.n_candidates.iloc[0]),
            'null_n_candidates_min': int(block.n_candidates.min()),
            'null_n_candidates_max': int(block.n_candidates.max())})
    result = pd.DataFrame(result)
    result['holm_adjusted_p'] = multipletests(result.exact_upper_tail_p, method='holm')[1]
    result['n_tests_in_family'] = len(result)
    return result


def write(name, table):
    table.to_csv(OUT / ('null_' + name), sep='\t', index=False)


def main():
    plan = ROOT / 'repro_null_plan.md'
    check('analysis_plan_exists_before_run', plan.exists())
    register(plan)
    lr = pd.read_csv(register(ROOT / 'literature/mouseconsensus.csv')).drop_duplicates(
        ['source_genesymbol', 'target_genesymbol'])
    all_metrics, all_gaps, all_dynamic, all_samples, all_effects, all_fixed = [], [], [], [], [], []
    all_allocations, coverage = [], []
    for config in CONFIGS:
        print('Loading', config, flush=True)
        data, lmats, manifests, perms = {}, {}, {}, {}
        for acc in DATASETS:
            folder, m = manifest(acc, config)
            manifests[acc] = m
            data[acc] = all_custom(acc, config, folder, m, lr)
            lmats[acc] = liana(acc, config, m)
            perms[acc] = allocations(m)
        common = data[DATASETS[0]]['scores'].index.intersection(data[DATASETS[1]]['scores'].index)
        for acc in DATASETS:
            for metric in LMETRICS:
                d = lmats[acc][metric]
                complete = d.index[np.isfinite(d.to_numpy()).all(axis=1)]
                common = common.intersection(complete)
        common = common.sort_values()
        check(config + '/nonempty_fixed_universe', len(common) > 0)
        fixed = common.to_frame(index=False); fixed.insert(0, 'config', config)
        all_fixed.append(fixed)
        scores = {}
        for acc in DATASETS:
            scores[acc] = {'custom_score': data[acc]['scores'].loc[common]}
            scores[acc].update({metric: lmats[acc][metric].loc[common] for metric in LMETRICS})
            for metric in METRICS:
                mat = scores[acc][metric]
                long = mat.rename_axis(columns='sample').stack().rename('value').reset_index()
                long.insert(0, 'metric', metric); long.insert(0, 'dataset', acc); long.insert(0, 'config', config)
                long = long.merge(manifests[acc], on='sample', validate='many_to_one')
                if metric in EXPR:
                    ranks = mat.rank(axis=0, method='average', pct=True).rename_axis(columns='sample').stack().rename(
                        'within_fixed_universe_percentile').reset_index()
                    long = long.merge(ranks, on=KEY + ['sample'], validate='one_to_one')
                else:
                    long['within_fixed_universe_percentile'] = np.nan
                all_samples.append(long)
                observed_effect = delta(mat, (manifests[acc].condition == 'MCAO').to_numpy())
                d = common.to_frame(index=False)
                d.insert(0, 'metric', metric); d.insert(0, 'dataset', acc); d.insert(0, 'config', config)
                d['observed_mcao_minus_sham'] = observed_effect
                all_effects.append(d)
            coverage.append({'config': config, 'dataset': acc, 'n_fixed_common_candidates': len(common),
                             'n_all_resource_gene_covered_edges': len(data[acc]['scores'])})
        # Precompute each cohort's effects/gates for all its allocations.
        fixed_effects, dynamic_effects, dynamic_gates = {}, {}, {}
        for acc in DATASETS:
            for p in perms[acc]:
                fixed_effects[acc, p['index']] = {metric: delta(scores[acc][metric], p['mask']) for metric in METRICS}
                dynamic_effects[acc, p['index']] = delta(data[acc]['scores'], p['mask'])
                dynamic_gates[acc, p['index']] = group_gate(data[acc], p['mask'])
        dacc, vacc = DATASETS
        for pid, (dp, vp) in enumerate(itertools.product(perms[dacc], perms[vacc])):
            meta = {'config': config, 'allocation_id': pid, 'discovery_allocation': dp['index'],
                    'validation_allocation': vp['index'], 'is_observed': dp['is_observed'] and vp['is_observed']}
            all_allocations.append(meta | {'discovery_mcao_animals': dp['mcao_animals'],
                                           'validation_mcao_animals': vp['mcao_animals']})
            results = {}
            for metric in METRICS:
                result = comparison(fixed_effects[dacc, dp['index']][metric], fixed_effects[vacc, vp['index']][metric])
                all_metrics.append(meta | {'metric': metric} | result)
                results[metric] = result
            for metric in EXPR:
                for stat in STATISTICS:
                    all_gaps.append(meta | {'expression_metric': metric, 'ranking_metric': 'magnitude_priority',
                                           'statistic': stat, 'value': results[metric][stat] - results['magnitude_priority'][stat],
                                           'n_candidates': len(common)})
            # Recompute disease-group detection eligibility in each permuted cohort.
            di = data[dacc]['scores'].index[dynamic_gates[dacc, dp['index']]]
            vi = data[vacc]['scores'].index[dynamic_gates[vacc, vp['index']]]
            dynamic_common = di.intersection(vi).sort_values()
            dpos = data[dacc]['scores'].index.get_indexer(dynamic_common)
            vpos = data[vacc]['scores'].index.get_indexer(dynamic_common)
            result = comparison(dynamic_effects[dacc, dp['index']][dpos], dynamic_effects[vacc, vp['index']][vpos])
            all_dynamic.append(meta | {'metric': 'custom_score', 'gate': .1} | result)
            if meta['is_observed']:
                # Confirm the original reported common-edge gate and effect statistics.
                old = []
                for acc in DATASETS:
                    base = ROOT / ('results' if config == 'primary' else 'revision_analysis/results') / acc
                    x = read(base / (config + '__communication.tsv.gz'))
                    old.append(x[x.eligible_10].set_index(KEY).score_difference)
                old_common = old[0].index.intersection(old[1].index).sort_values()
                check(config + '/observed_dynamic_gate_reproduces_original_keys', old_common.equals(dynamic_common), len(old_common))
                old_result = comparison(old[0].loc[old_common].to_numpy(), old[1].loc[old_common].to_numpy())
                check(config + '/observed_dynamic_gate_reproduces_original_statistics',
                      all(abs(old_result[s] - result[s]) < TOL for s in STATISTICS), result)
        print('Finished', config, 'fixed candidate count', len(common), flush=True)
    metric_table = pd.DataFrame(all_metrics)
    gap_table = pd.DataFrame(all_gaps)
    dynamic_table = pd.DataFrame(all_dynamic)
    index_cols = ['config', 'allocation_id', 'discovery_allocation', 'validation_allocation',
                  'is_observed', 'metric', 'n_candidates']
    metric_long = metric_table.melt(id_vars=index_cols, value_vars=STATISTICS, var_name='statistic', value_name='value')
    dynamic_long = dynamic_table.melt(id_vars=index_cols, value_vars=STATISTICS, var_name='statistic', value_name='value')
    observed_summary = summarize(metric_long, ['config', 'metric', 'statistic'], 'A_fixed_common_metrics')
    gaps_summary = summarize(gap_table, ['config', 'expression_metric', 'ranking_metric', 'statistic'], 'B_expression_minus_ranking')
    dynamic_summary = summarize(dynamic_long, ['config', 'metric', 'statistic'], 'C_exploratory_dynamic_gate')
    check('planned_test_family_sizes', [len(observed_summary), len(gaps_summary), len(dynamic_summary)] == [30, 24, 6])
    check('exact_probabilities_are_on_200_allocation_grid', all(np.allclose(x.exact_upper_tail_p * 200,
          np.round(x.exact_upper_tail_p * 200)) for x in [observed_summary, gaps_summary, dynamic_summary]))
    # Only after all calculations and assertions succeed do final outputs begin.
    OUT.mkdir(parents=True, exist_ok=True)
    write('fixed_common_summary.tsv', observed_summary)
    write('family_gap_summary.tsv', gaps_summary)
    write('dynamic_gate_summary.tsv', dynamic_summary)
    write('fixed_common_all_allocations.tsv.gz', metric_table)
    write('family_gap_all_allocations.tsv.gz', gap_table)
    write('dynamic_gate_all_allocations.tsv.gz', dynamic_table)
    write('animal_label_allocations.tsv', pd.DataFrame(all_allocations))
    write('fixed_candidate_keys.tsv', pd.concat(all_fixed, ignore_index=True))
    write('fixed_common_animal_scores.tsv.gz', pd.concat(all_samples, ignore_index=True))
    write('fixed_common_observed_effects.tsv.gz', pd.concat(all_effects, ignore_index=True))
    write('coverage.tsv', pd.DataFrame(coverage))
    metadata = {'status': 'complete', 'n_checks': len(CHECKS), 'checks': CHECKS,
                'source_sha256': INPUTS, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'n_joint_allocations': 200, 'observed_allocation_included': True,
                'randomization_unit': 'entire animal-level LR vector within each study',
                'zero_effect_absolute_tolerance': TOL, 'missing_scores_imputed': False,
                'fixed_universe_disease_label_dependent': False,
                'dynamic_gate_recomputed_per_allocation': True,
                'tail': 'one-sided upper', 'null_interval_is_effect_confidence_interval': False,
                'gap_null_tests_algorithm_equivalence': False,
                'multiplicity': {'A_fixed_common_metrics': 30, 'B_expression_minus_ranking': 24,
                                 'C_exploratory_dynamic_gate': 6},
                'software': {'numpy': np.__version__, 'pandas': pd.__version__}}
    (OUT / 'null_audit.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    (OUT / 'null_methods.md').write_text(plan.read_text(encoding='utf-8'), encoding='utf-8')
    print(observed_summary[observed_summary.metric.isin(['custom_score', 'magnitude_priority'])].to_string(index=False))
    print(dynamic_summary.to_string(index=False))


if __name__ == '__main__':
    main()
