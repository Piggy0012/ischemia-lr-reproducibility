"""Auditable animal-level summaries of the computational revision.

This script never fills missing LR scores with zero. The LIANA inputs must be
the final return_all_lrs=False exports, with full cell-type context upstream.
Run from any directory. A discovery-only self-check can use --datasets
GSE174574 --allow-partial --output-dir work/revision_summary_selfcheck.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
DATASETS = ['GSE174574', 'GSE245386']
CONFIGS = ['primary', 'singlet', 'reference_singlet', 'reference_only']
LIANA_CONFIGS = ['primary', 'reference_singlet']
KEY = ['sender', 'receiver', 'ligand', 'receptor']
RENAME = {'source': 'sender', 'target': 'receiver',
          'ligand_complex': 'ligand', 'receptor_complex': 'receptor'}
METRICS = {
    'lr_means': 'CellPhoneDB: ligand/receptor mean expression',
    'expr_prod': 'Connectome and NATMI: shared expression-product magnitude',
    'lr_logfc': 'log2FC: within-library cell-type-versus-rest enrichment',
    'spec_weight': 'NATMI: specificity weight',
    'lrscore': 'SingleCellSignalR: regularized expression score',
    'magnitude_rank': 'LIANA magnitude consensus, transformed to 1-rank',
}
EXAMPLES = [
    ('Astrocyte', 'Endothelial', 'Timp3', 'Kdr'),
    ('Endothelial', 'Astrocyte', 'Ptn', 'Ptprz1'),
    ('Endothelial', 'Astrocyte', 'Plat', 'Lrp1'),
    ('Astrocyte', 'Endothelial', 'Spp1', 'Itga5_Itgb1'),
    ('Astrocyte', 'Endothelial', 'Col4a1', 'Itga3_Itgb1'),
]


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path):
    return pd.read_csv(path, sep='\t')


def strict_bool(series):
    if series.dtype == bool:
        return series
    mapped = series.astype(str).str.lower().map({'true': True, 'false': False,
                                               '1': True, '0': False})
    if mapped.isna().any():
        raise ValueError(f'Invalid boolean values in {series.name}')
    return mapped.astype(bool)


def comparison(a, b):
    """Describe finite pairs; zeros have no positive/negative direction."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    direction = (a != 0) & (b != 0)
    same = ((a[direction] > 0) == (b[direction] > 0))
    rho = float(stats.spearmanr(a, b).statistic) if len(a) > 1 and np.ptp(a) and np.ptp(b) else np.nan
    return {'n_common': len(a), 'n_direction_evaluable': int(direction.sum()),
            'n_zero_effect': int((~direction).sum()), 'n_same_direction': int(same.sum()),
            'same_direction_fraction': float(same.mean()) if len(same) else np.nan,
            'spearman_rho': rho}


def effect_table(scores, conditions):
    """Scores are LR rows x all animal columns; incomplete rows are excluded."""
    case = np.asarray(conditions) == 'MCAO'
    if case.sum() < 2 or (~case).sum() < 2:
        raise ValueError('At least two animal libraries per condition are required')
    values = scores.to_numpy(float)
    complete = np.isfinite(values).all(axis=1)
    x = values[complete]
    result = scores.index.to_frame(index=False).loc[complete].reset_index(drop=True)
    result['n_sham'] = int((~case).sum())
    result['n_mcao'] = int(case.sum())
    result['mean_sham'] = x[:, ~case].mean(axis=1)
    result['mean_mcao'] = x[:, case].mean(axis=1)
    result['score_difference'] = result.mean_mcao - result.mean_sham
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        welch = stats.ttest_ind(x[:, case], x[:, ~case], axis=1, equal_var=False)
    result['animal_welch_p'] = welch.pvalue
    result['animal_welch_fdr'] = np.nan
    finite_p = np.isfinite(result.animal_welch_p)
    if finite_p.any():
        result.loc[finite_p, 'animal_welch_fdr'] = multipletests(
            result.loc[finite_p, 'animal_welch_p'], method='fdr_bh')[1]
    exceed = np.zeros(len(x), dtype=int)
    choices = list(itertools.combinations(range(len(case)), int(case.sum())))
    for selected in choices:
        mask = np.zeros(len(case), bool)
        mask[list(selected)] = True
        perm = x[:, mask].mean(axis=1) - x[:, ~mask].mean(axis=1)
        exceed += np.abs(perm) >= np.abs(result.score_difference.to_numpy()) - 1e-12
    result['animal_exact_permutation_p'] = exceed / len(choices)
    result['n_label_permutations'] = len(choices)
    result['n_bh_tests'] = int(finite_p.sum())
    result['status'] = np.where(finite_p, 'complete_case_tested', 'complete_case_welch_undefined')
    if len(result):
        assert np.all(result.animal_exact_permutation_p.between(1 / len(choices), 1))
        assert np.allclose(result.animal_exact_permutation_p * len(choices),
                           np.round(result.animal_exact_permutation_p * len(choices)))
    return result, complete


def json_safe(obj):
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def run(datasets, output, allow_partial=False):
    manifests, missing, required, inputs = {}, [], [], []
    for acc in datasets:
        mf = ROOT / 'processed' / acc / 'sample_cell_counts.tsv'
        m = read_table(mf)[['sample', 'condition']].drop_duplicates()
        assert not m['sample'].duplicated().any()
        assert set(m.condition) == {'Sham', 'MCAO'}
        manifests[acc] = m.sort_values('sample').reset_index(drop=True)
        required.append(mf)
        for config in CONFIGS:
            base = ROOT / ('results' if config == 'primary' else 'revision_analysis/results')
            required.append(base / acc / f'{config}__communication.tsv.gz')
        for config in LIANA_CONFIGS:
            for sample in m['sample']:
                stem = ROOT / 'revision_results/liana' / acc / f'{config}__{sample}'
                required.extend([Path(str(stem) + '.tsv.gz'), Path(str(stem) + '.json')])
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing and not allow_partial:
        raise FileNotFoundError('Inputs still pending; no final outputs written:\n' + '\n'.join(missing))
    for p in required:
        if p.exists():
            inputs.append({'path': str(p.relative_to(ROOT)), 'sha256': file_hash(p)})
    output.mkdir(parents=True, exist_ok=True)
    tables = {}

    # Full custom tables retain the original cohort-level expression gate.
    custom = {}
    example_rows = []
    for acc in datasets:
        for config in CONFIGS:
            base = ROOT / ('results' if config == 'primary' else 'revision_analysis/results')
            path = base / acc / f'{config}__communication.tsv.gz'
            if path.exists():
                d = read_table(path)
                assert not d.duplicated(KEY).any()
                for gate in ['05', '10', '20']:
                    d[f'eligible_{gate}'] = strict_bool(d[f'eligible_{gate}'])
                custom[acc, config] = d
            else:
                d = pd.DataFrame(columns=KEY)
            idx = d.set_index(KEY)
            for example in EXAMPLES:
                row = dict(zip(KEY, example)) | {'dataset': acc, 'config': config}
                if not path.exists():
                    row['status'] = 'input_pending'
                elif example not in idx.index:
                    row['status'] = 'absent_from_5pct_screen_or_resource_coverage'
                    row.update({'eligible_05': False, 'eligible_10': False, 'eligible_20': False})
                else:
                    row.update(idx.loc[example].to_dict())
                    row['status'] = 'eligible_10pct' if row['eligible_10'] else 'fails_10pct_gate'
                example_rows.append(row)
    tables['custom_fixed_examples'] = pd.DataFrame(example_rows)

    cross_rows, custom_details, within_rows = [], [], []
    for config in CONFIGS:
        if all((a, config) in custom for a in DATASETS):
            for gate in ['05', '10', '20']:
                d, v = [custom[a, config] for a in DATASETS]
                joined = d[d[f'eligible_{gate}']].merge(v[v[f'eligible_{gate}']], on=KEY,
                                                       suffixes=('_discovery', '_validation'))
                cross_rows.append({'config': config, 'expression_gate': float(gate) / 100,
                                   'n_discovery_eligible': int(d[f'eligible_{gate}'].sum()),
                                   'n_validation_eligible': int(v[f'eligible_{gate}'].sum()),
                                   **comparison(joined.score_difference_discovery, joined.score_difference_validation)})
                joined['same_direction'] = (joined.score_difference_discovery * joined.score_difference_validation) > 0
                joined.insert(0, 'expression_gate', float(gate) / 100)
                joined.insert(0, 'selection', config)
                custom_details.append(joined)
    for (acc, config), d in custom.items():
        if config != 'primary' and (acc, 'primary') in custom:
            base = custom[acc, 'primary']
            joined = base[base.eligible_10].merge(d[d.eligible_10], on=KEY,
                                                suffixes=('_original', '_revised'))
            within_rows.append({'dataset': acc, 'config': config, 'comparator': 'primary',
                                **comparison(joined.score_difference_original, joined.score_difference_revised)})
    tables['custom_cross_cohort_summary'] = pd.DataFrame(cross_rows)
    tables['custom_cross_cohort_edges'] = pd.concat(custom_details, ignore_index=True) if custom_details else pd.DataFrame(columns=['selection', 'expression_gate'] + KEY)
    tables['custom_vs_original_within_cohort'] = pd.DataFrame(within_rows)

    # LIANA is evaluated separately in each animal, including complete-case status.
    liana_effects, liana_samples, eligibility_rows, liana_examples = [], [], [], []
    liana_audits, contrast_meta, effect_lookup = [], [], {}
    for acc in datasets:
        m = manifests[acc]
        sample_order = m['sample'].tolist()
        for config in LIANA_CONFIGS:
            per_sample, present = [], set()
            for sample, condition in m.itertuples(index=False, name=None):
                stem = ROOT / 'revision_results/liana' / acc / f'{config}__{sample}'
                path, ap = Path(str(stem) + '.tsv.gz'), Path(str(stem) + '.json')
                if not path.exists() or not ap.exists():
                    for example in EXAMPLES:
                        eligibility_rows.append(dict(zip(KEY, example)) | {'dataset': acc, 'config': config,
                            'sample': sample, 'condition': condition, 'status': 'input_pending'})
                    continue
                audit = json.loads(ap.read_text(encoding='utf-8'))
                assert audit['parameters']['return_all_lrs'] is False, 'Stale ineligible-filled LIANA export'
                assert audit['computed_all_cell_type_pairs'] is True, 'Full context network required'
                d = read_table(path).rename(columns=RENAME)
                assert not d.duplicated(KEY).any()
                assert set(d['sample']) <= {sample} and set(d.condition) <= {condition}
                assert set(d.dataset) <= {acc} and set(d.config) <= {config}
                assert strict_bool(d.expression_eligible).all()
                assert set(METRICS) <= set(d.columns)
                assert d.magnitude_rank.dropna().between(0, 1).all()
                d['magnitude_priority'] = 1 - d.magnitude_rank
                present.add(sample)
                per_sample.append(d)
                liana_audits.append({'dataset': acc, 'config': config, 'sample': sample,
                    'condition': condition, 'n_exported_eligible_rows': len(d),
                    'all_pair_rows': audit['all_pair_rows'], 'n_context_cells': audit['n_context_cells'],
                    'n_context_types': len(audit['context_counts']), 'liana_version': audit['liana_version']})
                idx = d.set_index(KEY)
                for example in EXAMPLES:
                    rec = dict(zip(KEY, example)) | {'dataset': acc, 'config': config,
                                                   'sample': sample, 'condition': condition}
                    if example in idx.index:
                        rec.update(idx.loc[example].to_dict())
                        rec.update({'status': 'reported_expression_eligible', 'reported_eligible': True})
                    else:
                        rec.update({'status': 'not_reported_by_expression_or_resource_filter', 'reported_eligible': False})
                    eligibility_rows.append(rec)
            if not per_sample:
                continue
            observed = pd.concat(per_sample, ignore_index=True)
            liana_samples.append(observed)
            if present != set(sample_order):
                contrast_meta.append({'dataset': acc, 'config': config, 'status': 'input_pending',
                                      'missing_samples': sorted(set(sample_order) - present)})
                continue
            for metric in METRICS:
                value_col = 'magnitude_priority' if metric == 'magnitude_rank' else metric
                scores = observed.pivot(index=KEY, columns='sample', values=value_col).reindex(columns=sample_order)
                # No fillna: absent entries are not measured zero expression.
                effect, complete = effect_table(scores, m.condition.to_numpy())
                effect.insert(0, 'metric', metric)
                effect.insert(0, 'config', config)
                effect.insert(0, 'dataset', acc)
                liana_effects.append(effect)
                effect_lookup[acc, config, metric] = effect
                contrast_meta.append({'dataset': acc, 'config': config, 'metric': metric,
                    'status': 'complete', 'n_union_reported_edges': len(scores),
                    'n_complete_case_edges': int(complete.sum()), 'n_incomplete_edges': int((~complete).sum()),
                    'n_bh_tests': int(effect.animal_welch_p.notna().sum()),
                    'n_welch_fdr_05': int((effect.animal_welch_fdr < .05).sum()),
                    'min_exact_p': effect.animal_exact_permutation_p.min() if len(effect) else np.nan})
                idx = effect.set_index(KEY)
                for example in EXAMPLES:
                    rec = dict(zip(KEY, example)) | {'dataset': acc, 'config': config, 'metric': metric}
                    values = scores.loc[example] if example in scores.index else pd.Series(np.nan, index=sample_order)
                    finite_values = np.isfinite(values).to_numpy()
                    case = m.condition.to_numpy() == 'MCAO'
                    rec.update({'n_expected_animals': len(sample_order),
                                'n_observed_finite': int(finite_values.sum()),
                                'n_observed_sham': int((finite_values & ~case).sum()),
                                'n_observed_mcao': int((finite_values & case).sum()),
                                'missing_or_nonfinite_samples': '|'.join(values.index[~finite_values])})
                    if example in idx.index:
                        rec.update(idx.loc[example].to_dict())
                    else:
                        rec['status'] = 'not_complete_case'
                    liana_examples.append(rec)
    tables['liana_per_animal_observed_scores'] = pd.concat(liana_samples, ignore_index=True) if liana_samples else pd.DataFrame(columns=KEY)
    tables['liana_per_animal_audit'] = pd.DataFrame(liana_audits)
    tables['liana_fixed_examples_per_animal_eligibility'] = pd.DataFrame(eligibility_rows)
    tables['liana_animal_level_effects'] = pd.concat(liana_effects, ignore_index=True) if liana_effects else pd.DataFrame(columns=KEY)
    tables['liana_contrast_coverage'] = pd.DataFrame(contrast_meta)
    tables['liana_fixed_examples_animal_effects'] = pd.DataFrame(liana_examples)

    liana_cross, liana_detail, method_custom, method_custom_detail = [], [], [], []
    for config in LIANA_CONFIGS:
        for metric in METRICS:
            if all((a, config, metric) in effect_lookup for a in DATASETS):
                d, v = [effect_lookup[a, config, metric] for a in DATASETS]
                joined = d.merge(v, on=KEY, suffixes=('_discovery', '_validation'))
                liana_cross.append({'config': config, 'metric': metric,
                    'n_discovery_complete_edges': len(d), 'n_validation_complete_edges': len(v),
                    **comparison(joined.score_difference_discovery, joined.score_difference_validation)})
                joined['same_direction'] = (joined.score_difference_discovery * joined.score_difference_validation) > 0
                joined.insert(0, 'selection', config)
                joined.insert(0, 'score_metric', metric)
                liana_detail.append(joined)
            for acc in datasets:
                if (acc, config, metric) not in effect_lookup or (acc, config) not in custom:
                    continue
                c = custom[acc, config]
                joined = effect_lookup[acc, config, metric].merge(c[c.eligible_10], on=KEY,
                                                                 suffixes=('_liana', '_custom'))
                method_custom.append({'dataset': acc, 'config': config, 'metric': metric,
                    'custom_gate': .1, **comparison(joined.score_difference_liana, joined.score_difference_custom)})
                joined['same_direction'] = (joined.score_difference_liana * joined.score_difference_custom) > 0
                method_custom_detail.append(joined)
    tables['liana_cross_cohort_summary'] = pd.DataFrame(liana_cross)
    tables['liana_cross_cohort_edges'] = pd.concat(liana_detail, ignore_index=True) if liana_detail else pd.DataFrame(columns=KEY)
    tables['liana_vs_custom_direction_summary'] = pd.DataFrame(method_custom)
    tables['liana_vs_custom_direction_edges'] = pd.concat(method_custom_detail, ignore_index=True) if method_custom_detail else pd.DataFrame(columns=KEY)

    for name, table in tables.items():
        # Compressed edges retain every evaluated row; short summaries remain readable.
        suffix = '.tsv.gz' if name.endswith('_edges') or name in ['liana_per_animal_observed_scores', 'liana_animal_level_effects'] else '.tsv'
        table.to_csv(output / (name + suffix), sep='\t', index=False)
    complete = not missing and set(datasets) == set(DATASETS)
    summary = {'status': 'complete' if complete else 'partial_not_for_final_manuscript',
        'datasets_requested': datasets, 'pending_inputs': missing,
        'metrics': METRICS, 'fixed_examples': [dict(zip(KEY, e)) for e in EXAMPLES],
        'custom_cross_cohort': cross_rows, 'custom_vs_original_within_cohort': within_rows,
        'liana_contrast_coverage': contrast_meta, 'liana_cross_cohort': liana_cross,
        'liana_vs_custom': method_custom, 'table_rows': {n: len(t) for n, t in tables.items()},
        'script_sha256': file_hash(Path(__file__)), 'input_files': inputs,
        'inference_unit': 'animal library', 'missing_scores_imputed': False,
        'liana_internal_pvalues_used_for_disease_inference': False}
    (output / 'summary.json').write_text(json.dumps(json_safe(summary), ensure_ascii=False,
                                                   indent=2, allow_nan=False), encoding='utf-8')
    (output / 'methodsREADME.md').write_text(METHODS_README, encoding='utf-8')
    print(json.dumps({'status': summary['status'], 'output': str(output),
                      'table_rows': summary['table_rows']}, ensure_ascii=False), flush=True)
    return summary


METHODS_README = '''# 修订结果表的统计口径

这些表为同一公开数据的计算敏感性和方法复核；算法一致不是新增独立队列。

custom 表沿用每研究原先保存的5%、10%、20%表达门槛，疾病效应为MCAO减Sham的动物级共可用性均值差。跨队列只比较双方均满足相应门槛的候选；另列三种新选择与原primary在各研究内的共同10%候选一致性。固定四个原示例及Col4a1–Itga3_Itgb1全部保留。源custom文件未输出未过5%筛选或资源不覆盖的候选，故此类的效应为空，不补零。Welch FDR直接沿用对应门槛全候选集合的结果。

LIANA仅接受最终return_all_lrs=False、每文库全细胞类型网络计算后导出的双向目标结果。原始每动物全方法数值单独保存，cellphone_pvals、specificity_rank不是动物间疾病差异P值。每研究、细胞选择、指标分别以所有动物都有有限已报告分数的候选构成complete-case集合；表达筛选或资源筛选导致的未报告和任何非有限值均不补零。固定示例表逐动物列出资格与缺失，可据此核查筛选造成的分母变化。未报告不等同于生物学上没有通讯。complete-case筛选依赖观测表达，因此其检验是条件于可评估集合的探索性推断。

疾病比较指标为lr_means、expr_prod、lr_logfc、spec_weight、lrscore及1-magnitude_rank。前五列覆盖CellPhoneDB、Connectome、log2FC、NATMI、SingleCellSignalR五个实际方法；expr_prod由Connectome和NATMI共用，不能算两份独立证据。lr_logfc本身是文库内细胞类型相对其它类型的表达富集；其MCAO-Sham差是富集指标差，不是疾病基因log2FC。1-magnitude_rank使数值增大表示相对优先级提高；各文库全网络候选空间可不同，它不是绝对通讯强度，亦不可把两个rank当相同标尺的实验测量。

每个LIANA疾病效应是MCAO减Sham的动物均值差。Welch检验按每研究×选择×指标的双方向完整候选集合做BH校正；没有对6个相关指标做一个共同校正，不能从中择取最小P值作确认性结论。恒定或退化数据若产生未定义Welch P值则保留NA，不写成0或显著；BH仅作用于有限P值并报告n_bh_tests。另穷举整只动物的处理标签，以均值差绝对值为双侧统计量，计算(|置换差| >= |观察差|-1e-12)的比例；不是细胞标签置换。3对3共20种、3对2共10种标签分配，离散分辨率限制保持不变。精确置换基于可交换性假设，仍不消除队列设计限制。

一致性表中的n_common为共同有限效应数；n_direction_evaluable排除任一效应恰好等于0的条目，同向比例以此为分母；Spearman使用所有共同有限效应。共享配体或受体使边之间不独立，因此不输出基于独立边假定的一致率P值或相关显著性。LIANA与custom只在共同可评估候选上作方向描述。低表达导致完整病例候选集合收缩必须与一致率同时报告。

summary.json保存输入SHA-256、脚本SHA-256、每表行数、各分母与待生成文件。partial_not_for_final_manuscript只用于过程自检；正式运行须全部输入存在，默认缺任何必需文件即终止且不写最终表。
'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets', nargs='+', choices=DATASETS, default=DATASETS)
    parser.add_argument('--output-dir', type=Path, default=ROOT.parent / 'outputs/revised/tables')
    parser.add_argument('--allow-partial', action='store_true')
    args = parser.parse_args()
    run(args.datasets, args.output_dir.resolve(), args.allow_partial)
