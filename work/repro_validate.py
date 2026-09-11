"""Read-only scientific validation of cached v2 outputs; no cell matrices or fits.

python work/repro_validate.py [--require-complete] [--report-dir DIRECTORY]
Exit 0: completed stages pass (optional stages may be explicitly pending).
Exit 1: a completed-stage check failed. Exit 2: --require-complete has pending stages.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import re
import sys
import traceback

import numpy as np
import pandas as pd
from scipy import stats

WORK = Path(__file__).resolve().parent
BASE = WORK.parent
TABLES = BASE / 'outputs/reproducibility_v2/tables'
KEY = ['sender', 'receiver', 'ligand', 'receptor']
CONFIGS = {'primary': 187, 'reference_singlet': 174}
GROUPS = {'GSE174574': {'MCAO': 3, 'Sham': 3}, 'GSE245386': {'MCAO': 2, 'Sham': 3}}
METRICS = ['custom_score', 'lr_means', 'expr_prod', 'lrscore', 'magnitude_priority']
STATS = ['same_direction_all_fraction', 'same_direction_nonzero_fraction', 'spearman_rho']
TOL = 1e-12
REQUIRED = {
    'null': ['audit.json', 'animal_label_allocations.tsv', 'coverage.tsv', 'fixed_candidate_keys.tsv',
             'fixed_common_animal_scores.tsv.gz', 'fixed_common_observed_effects.tsv.gz',
             'fixed_common_all_allocations.tsv.gz', 'fixed_common_summary.tsv',
             'family_gap_all_allocations.tsv.gz', 'family_gap_summary.tsv',
             'dynamic_gate_all_allocations.tsv.gz', 'dynamic_gate_summary.tsv'],
    'null_conditional': ['audit.json', 'all_allocations.tsv.gz', 'summary.tsv', 'methods.md', 'explanation_zh.md'],
    'rank': ['audit.json', 'comparison_summary.tsv', 'figure_source_comparison.tsv', 'observed_effects.tsv.gz',
             'per_animal_representations.tsv.gz', 'network_context_per_animal.tsv', 'network_context_cell_counts.tsv',
             'native_target_ties_per_animal.tsv', 'native_zero_effects.tsv', 'network_context_ranges.tsv',
             'implementation_four_row_probe.tsv'],
    'rank_diagnostic_raw': ['audit.json', 'summary.tsv', 'animal_scores.tsv.gz', 'effects.tsv.gz', 'ties.tsv'],
    'rank_background': ['audit.json', 'summary.tsv', 'animal_scores.tsv.gz', 'effects.tsv.gz',
                        'global_keys.tsv.gz', 'monotonic_redundancy.tsv', 'context.tsv', 'zero_effects.tsv'],
    'rank_upstream': ['audit.json', 'per_library_equality.tsv', 'animal_scores.tsv.gz', 'effects.tsv.gz', 'summary.tsv'],
    'sorted': ['descriptive_audit.json', 'all_gene_descriptive_effects.tsv.gz', 'fixed_component_effects_ci.tsv'],
    'third': ['coavailability_audit.json', 'all_covered_coavailability_effects.tsv.gz',
              'all_covered_coavailability_library_scores.tsv.gz', 'cell_type_library_counts.tsv',
              'cross_cohort_concordance_edges.tsv.gz', 'cross_cohort_concordance_summary.tsv',
              'fixed_examples_component_expression.tsv', 'fixed_examples_effects_eligibility.tsv',
              'fixed_examples_library_scores.tsv', 'gene_feature_coverage.tsv.gz', 'library_metadata.tsv'],
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def resolve_recorded_path(value):
    """Relocate original Windows audit paths without depending on the old checkout."""
    value = str(value).replace('\\', '/')
    if '/site-packages/liana/' in value:
        spec = importlib.util.find_spec('liana')
        if spec is None:
            raise FileNotFoundError('Install the archived LIANA dependency to verify its source hash')
        return Path(spec.origin).parent / value.split('/site-packages/liana/', 1)[1]
    # A relocated release can itself sit below an older work/ directory. Use
    # the innermost archive-root marker, not a fixed preference for work/.
    markers = [(value.rfind('/work/'), '/work/', WORK),
               (value.rfind('/outputs/'), '/outputs/', BASE / 'outputs')]
    index, marker, root = max(markers, key=lambda item: item[0])
    if index >= 0:
        return root / value[index + len(marker):]
    if value.startswith(('work/', 'outputs/')):
        return BASE / value
    if (WORK / value).exists():
        return WORK / value
    if (TABLES / value).exists():
        return TABLES / value
    return WORK / value


def holm(pvalues):
    p = np.asarray(pvalues, float)
    order = np.argsort(p, kind='stable')
    corrected = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]))
    result = np.empty_like(corrected)
    result[order] = corrected
    return result


class Validator:
    def __init__(self):
        self.checks = []
        self.warnings = []
        self.pending = []
        self.stages = []
        self.files = {}
        self._tables = {}

    def check(self, name, condition, detail=None):
        ok = bool(condition)
        self.checks.append({'check': name, 'passed': ok, 'detail': detail})
        if not ok:
            raise AssertionError(name + ': ' + str(detail))

    def read(self, filename):
        if filename not in self._tables:
            path = TABLES / filename
            self.check('exists/' + filename, path.is_file() and path.stat().st_size > 0)
            # Preserve the exact serialized binary float. The default C parser
            # can split/merge almost-equal percentile effects and change rho.
            self._tables[filename] = pd.read_csv(path, sep='\t', float_precision='round_trip')
        return self._tables[filename]

    def audit(self, filename, script=None):
        path = TABLES / filename
        a = json.loads(path.read_text(encoding='utf-8'))
        if 'status' in a:
            self.check(filename + '/complete', a['status'].startswith('complete'))
        if 'checks' in a:
            self.check(filename + '/recorded_checks', all(x['passed'] for x in a['checks']))
            self.check(filename + '/recorded_check_count', a['n_checks'] == len(a['checks']))
        for field in ['source_sha256', 'input_sha256']:
            if isinstance(a.get(field), dict):
                for recorded, expected in a[field].items():
                    actual = resolve_recorded_path(recorded)
                    self.check(filename + '/input_hash/' + recorded, actual.exists() and sha(actual) == expected)
        if script and a.get('script_sha256'):
            self.check(filename + '/script_hash', sha(WORK / script) == a['script_sha256'])
        return a

    def stage(self, name, fn):
        start = len(self.checks)
        try:
            fn()
            status = 'passed'
        except Exception as exc:
            if not self.checks or self.checks[-1]['passed']:
                self.checks.append({'check': name + '/exception', 'passed': False,
                                    'detail': f'{type(exc).__name__}: {exc}'})
            status = 'failed'
        self.stages.append({'stage': name, 'status': status, 'n_checks': len(self.checks) - start})

    def inventory(self):
        for prefix, suffixes in REQUIRED.items():
            for suffix in suffixes:
                path = TABLES / f'{prefix}_{suffix}'
                self.check('required/' + path.name, path.is_file() and path.stat().st_size > 0)
                self.files[str(path.relative_to(BASE))] = sha(path)
                if not path.name.endswith(('.tsv', '.tsv.gz')):
                    continue
                d = self.read(path.name)
                numeric = d.select_dtypes(include=np.number)
                self.check(path.name + '/no_infinities', not np.isinf(numeric.to_numpy()).any())
                allowed = {}
                if path.name == 'null_fixed_common_animal_scores.tsv.gz':
                    allowed['within_fixed_universe_percentile'] = d.metric.eq('magnitude_priority')
                elif path.name == 'null_conditional_all_allocations.tsv.gz':
                    allowed['gate'] = d.analysis.eq('fixed_common')
                elif path.name == 'third_library_metadata.tsv':
                    allowed['n_independent_animals'] = np.ones(len(d), dtype=bool)
                for col in d:
                    self.check(path.name + '/missing/' + col,
                               not d[col].isna().any() or
                               (col in allowed and np.array_equal(d[col].isna().to_numpy(), np.asarray(allowed[col]))))
        self.check('historical_sorted_tests_retained', (WORK / 'sorted_rna/diff_expression.tsv.gz').is_file())

    def animal_effects(self, scores_name, effects_name, effect_column, extra=()):
        scores = self.read(scores_name)
        effects = self.read(effects_name)
        groups = ['config', 'dataset', 'metric', *extra]
        self.check(scores_name + '/unique_sample_edge_rows',
                   not scores.duplicated(groups + KEY + ['sample']).any())
        for group, block in scores.groupby(groups, sort=False):
            config, dataset = group[:2]
            m = block[['sample', 'condition']].drop_duplicates().sort_values('sample')
            self.check(str(group) + '/one_condition_per_animal', m['sample'].is_unique)
            self.check(str(group) + '/expected_animal_counts', m.condition.value_counts().to_dict() == GROUPS[dataset])
            matrix = block.pivot(index=KEY, columns='sample', values='value').reindex(columns=m['sample'])
            self.check(str(group) + '/common_finite_matrix', len(matrix) == CONFIGS[config] and np.isfinite(matrix).all().all())
            case = m.condition.to_numpy() == 'MCAO'
            values = matrix.to_numpy()
            delta = values[:, case].mean(axis=1) - values[:, ~case].mean(axis=1)
            delta[np.abs(delta) <= TOL] = 0
            target = effects
            for col, value in zip(groups, group):
                target = target[target[col] == value]
            target = target.set_index(KEY).reindex(matrix.index)
            self.check(str(group) + '/equal_animal_effect_reconstruction',
                       len(target) == len(delta) and np.allclose(target[effect_column], delta, atol=TOL, rtol=0))

    def summary_from_effects(self, effects_name, summary_name, effect_col, extra=()):
        e, s = self.read(effects_name), self.read(summary_name)
        groups = ['config', 'metric', *extra]
        for group, block in e.groupby(groups, sort=False):
            mat = block.pivot(index=KEY, columns='dataset', values=effect_col)
            self.check(str(group) + '/identical_cohort_candidates',
                       list(sorted(mat.columns)) == list(sorted(GROUPS)) and np.isfinite(mat).all().all())
            a, b = (mat[x].to_numpy(copy=True) for x in GROUPS)
            a[np.abs(a) <= TOL] = 0
            b[np.abs(b) <= TOL] = 0
            nz = (a != 0) & (b != 0)
            same = nz & (np.sign(a) == np.sign(b))
            row = s
            for col, value in zip(groups, group):
                row = row[row[col] == value]
            self.check(str(group) + '/one_summary', len(row) == 1)
            row = row.iloc[0]
            expected = {'n_candidates': len(a), 'n_nonzero_both': int(nz.sum()),
                        'n_same_direction': int(same.sum()), 'same_direction_all_fraction': same.mean(),
                        'same_direction_nonzero_fraction': same.sum() / nz.sum(),
                        'spearman_rho': stats.spearmanr(a, b).statistic}
            for col, value in expected.items():
                self.check(str(group) + '/' + col, np.isclose(row[col], value, atol=TOL, rtol=0))

    def nulls(self):
        a = self.audit('null_audit.json', 'repro_null_baseline.py')
        self.check('null/animal_unit', a['randomization_unit'] == 'entire animal-level LR vector within each study')
        self.check('null/design', a['n_joint_allocations'] == 200 and a['observed_allocation_included'] and
                   not a['missing_scores_imputed'] and not a['fixed_universe_disease_label_dependent'] and
                   a['dynamic_gate_recomputed_per_allocation'] and not a['gap_null_tests_algorithm_equivalence'])
        labels = self.read('null_animal_label_allocations.tsv')
        for config, block in labels.groupby('config'):
            self.check('allocations/' + config, len(block) == 200 and block.is_observed.sum() == 1 and
                       not block.duplicated(['discovery_allocation', 'validation_allocation']).any() and
                       block.discovery_allocation.nunique() == 20 and block.validation_allocation.nunique() == 10)
            self.check('allocations/' + config + '/case_counts',
                       block.discovery_mcao_animals.str.split('|', regex=False).map(len).eq(3).all() and
                       block.validation_mcao_animals.str.split('|', regex=False).map(len).eq(2).all())
        keys = self.read('null_fixed_candidate_keys.tsv')
        self.check('null/fixed_counts', keys.groupby('config').size().to_dict() == CONFIGS)
        self.check('null/fixed_keys_unique', not keys.duplicated(['config'] + KEY).any())
        self.animal_effects('null_fixed_common_animal_scores.tsv.gz', 'null_fixed_common_observed_effects.tsv.gz', 'observed_mcao_minus_sham')
        specs = [('fixed_common', 'A_fixed_common_metrics', 30), ('family_gap', 'B_expression_minus_ranking', 24),
                 ('dynamic_gate', 'C_exploratory_dynamic_gate', 6)]
        for prefix, family, n in specs:
            d = self.read(f'null_{prefix}_all_allocations.tsv.gz')
            s = self.read(f'null_{prefix}_summary.tsv')
            self.check(family + '/size_and_label', len(s) == n and s.test_family.eq(family).all() and s.n_tests_in_family.eq(n).all())
            for r in s.itertuples(index=False):
                b = d[d.config == r.config]
                if prefix == 'family_gap':
                    b = b[(b.expression_metric == r.expression_metric) & (b.ranking_metric == r.ranking_metric) & (b.statistic == r.statistic)]
                    values = b.value.to_numpy()
                else:
                    b = b[b.metric == r.metric]
                    values = b[r.statistic].to_numpy()
                self.check(family + '/full_distribution', len(b) == 200 and b.is_observed.sum() == 1)
                observed = values[b.is_observed.to_numpy()][0]
                k = int((values >= observed - TOL).sum())
                expected = [observed, values.mean(), np.median(values), *np.quantile(values, [.025, .975]), k / 200]
                actual = [r.observed, r.null_mean, r.null_median, r.null_q025, r.null_q975, r.exact_upper_tail_p]
                self.check(family + '/tail_and_reference', np.allclose(expected, actual, atol=TOL, rtol=0) and k == r.n_greater_or_equal_observed)
            self.check(family + '/independent_Holm', np.allclose(holm(s.exact_upper_tail_p), s.holm_adjusted_p, atol=TOL, rtol=0))
        self.check('null/exact_family_registry', a['multiplicity'] == {x[1]: x[2] for x in specs})

    def conditional(self):
        a = self.audit('null_conditional_audit.json', 'repro_conditional_null.py')
        self.check('conditional/exploratory_status', not a['new_confirmatory_test'] and not a['original_families_changed'] and
                   a['joint_results_known_before_diagnostic_plan'] and a['conditional_pvalues_are_unadjusted_diagnostics'])
        d, s = self.read('null_conditional_all_allocations.tsv.gz'), self.read('null_conditional_summary.tsv')
        self.check('conditional/72_rows', len(s) == 72 and not s.conditional_p_adjusted_for_multiple_testing.any())
        for r in s.itertuples(index=False):
            b = d[(d.analysis == r.analysis) & (d.config == r.config) & (d.metric == r.metric) & (d.randomized_cohort == r.randomized_cohort)]
            joint = self.read('null_' + r.analysis + '_all_allocations.tsv.gz')
            joint = joint[(joint.config == r.config) & (joint.metric == r.metric)]
            obs = joint[joint.is_observed].iloc[0]
            fixed = 'validation_allocation' if r.randomized_cohort == 'GSE174574' else 'discovery_allocation'
            expected_ids = set(joint[joint[fixed] == obs[fixed]].allocation_id)
            n = 20 if r.randomized_cohort == 'GSE174574' else 10
            self.check('conditional/exact_subset', len(b) == n and set(b.allocation_id) == expected_ids and b.is_observed.sum() == 1)
            values = b[r.statistic].to_numpy()
            p = float((values >= r.observed - TOL).mean())
            self.check('conditional/tail', np.isclose(p, r.conditional_unadjusted_upper_tail_p, atol=TOL, rtol=0) and r.minimum_attainable_upper_tail_p == 1 / n)

    def ranks(self):
        for filename, script in [('rank_audit.json', 'repro_rank_mechanism.py'), ('rank_background_audit.json', 'repro_rank_background.py')]:
            a = self.audit(filename, script)
            self.check(filename + '/descriptive', not a.get('pvalues_computed', a.get('new_pvalues', True)))
        for stem, extra in [('rank', ['representation']), ('rank_diagnostic_raw', []), ('rank_background', []), ('rank_upstream', [])]:
            scores = 'rank_per_animal_representations.tsv.gz' if stem == 'rank' else stem + '_animal_scores.tsv.gz'
            effects = 'rank_observed_effects.tsv.gz' if stem == 'rank' else stem + '_effects.tsv.gz'
            summary = 'rank_comparison_summary.tsv' if stem == 'rank' else stem + '_summary.tsv'
            self.animal_effects(scores, effects, 'mcao_minus_sham', extra)
            self.summary_from_effects(effects, summary, 'mcao_minus_sham', extra)
        keys = self.read('rank_background_global_keys.tsv.gz')
        self.check('background/global_counts', keys.groupby('config').size().to_dict() == {'primary': 7675, 'reference_singlet': 7109})
        targets = self.read('null_fixed_candidate_keys.tsv')
        self.check('background/includes_fixed_targets', len(targets.merge(keys, on=['config'] + KEY, validate='1:1')) == len(targets))
        z = self.read('rank_background_zero_effects.tsv')
        self.check('background/zero_effect_provenance', z.n_zero_disease_effect.eq(z.n_zero_effect_with_all_animals_at_priority_zero).all())
        u = self.audit('rank_upstream_audit.json', 'repro_upstream_correction_audit.py')
        self.check('upstream/attribution_scope', u['n_library_config_audits'] == 22 and
                   u['known_upstream_fix_not_a_new_bug_discovery'] and not u['new_pvalues'] and
                   not u['new_full_pipeline_run'] and not u['installed_package_modified'])
        self.check('upstream/exact_functions_and_outputs', u['function_bodies_AST_identical'] and
                   u['all_upstream_vs_original_unique_diagnostic_bitwise_equal'] and
                   u['maximum_absolute_difference_from_original_unique_diagnostic'] == 0)

    def diagnostic_networks(self, source='raw'):
        name = f'rank_diagnostic_{source}_audit.json'
        a = self.audit(name)
        self.check(name + '/full_22', a['n_library_config_audits'] == 22 and set(a['configs']) == set(CONFIGS))
        self.check(name + '/no_new_inference', not a['inferential_pvalues_computed'])
        expected_script = sha(WORK / 'repro_liana_rank_diagnostic_raw_run.py') if source == 'raw' else None
        found = set()
        for rec in a['audits']:
            identity = (rec['dataset'], rec['sample'], rec['config'])
            self.check(name + '/unique_library_config', identity not in found)
            found.add(identity)
            folder = WORK / 'repro_liana_rank_diagnostic' / source / rec['dataset']
            stem = folder / (rec['config'] + '__' + rec['sample'])
            local = json.loads(Path(str(stem) + '.json').read_text(encoding='utf-8'))
            self.check(str(identity) + '/audit_matches_summary', local == rec)
            self.check(str(identity) + '/complete', rec['status'] == 'complete' and rec['count_source'] == source and rec['computed_all_cell_type_pairs'])
            self.check(str(identity) + '/unaltered_package', not rec['package_modified'] and not rec['diagnostic_is_official_package_release'])
            self.check(str(identity) + '/native_reconstruction', rec['native_full_network_reconstruction_max_abs_error'] <= TOL)
            if source == 'raw':
                self.check(str(identity) + '/raw_run_snapshot_hash', rec['script_sha256'] == expected_script)
                self.check(str(identity) + '/raw_prior_exact', rec['raw_target_scores_verified_against_prior'] and
                           max(rec['raw_target_prior_max_abs_errors'].values()) <= TOL)
            for suffix, field in [('__full_network.tsv.gz', 'full_network_sha256'), ('__target.tsv.gz', 'target_sha256')]:
                path = Path(str(stem) + suffix)
                self.check(str(identity) + '/' + field, path.exists() and sha(path) == rec[field])
                d = pd.read_csv(path, sep='\t', usecols=['source', 'target', 'ligand_complex', 'receptor_complex',
                                                        'lr_means', 'expr_prod', 'lrscore', 'magnitude_rank',
                                                        'diagnostic_unique_column_magnitude_rank'])
                self.check(str(identity) + suffix + '/keys_finite', not d.duplicated(list(d.columns[:4])).any() and
                           np.isfinite(d.select_dtypes(include=np.number)).all().all())
        for config in CONFIGS:
            for acc, counts in GROUPS.items():
                self.check(name + '/' + config + '/' + acc, sum(x[0] == acc and x[2] == config for x in found) == sum(counts.values()))

    def sorted(self):
        a = self.audit('sorted_descriptive_audit.json', 'repro_sorted_descriptive.py')
        self.check('sorted/no_new_significance', not a['new_p_values_computed'] and not a['significance_classification'])
        self.check('sorted/source_hash', sha(WORK / 'sorted_rna/astro_endothelial_normalized.tsv.gz') == a['source_sha256'])
        d = self.read('sorted_all_gene_descriptive_effects.tsv.gz')
        self.check('sorted/units_and_intervals', len(d) == 49914 and d.n_paired_pools.eq(6).all() and
                   (d.ci_low <= d.paired_log2_effect).all() and (d.ci_high >= d.paired_log2_effect).all())
        old = pd.read_csv(WORK / 'sorted_rna/diff_expression.tsv.gz', sep='\t')
        j = old.merge(d, on=['cell_type', 'gene'], suffixes=('_old', '_new'), validate='1:1')
        self.check('sorted/historical_effects_preserved', len(j) == len(old) and
                   np.allclose(j.log2_difference, j.paired_log2_effect, atol=1e-12, rtol=0) and
                   np.allclose(j.ci_low_old, j.ci_low_new, atol=1e-12, rtol=0) and
                   np.allclose(j.ci_high_old, j.ci_high_new, atol=1e-12, rtol=0))

    def third(self):
        a = self.audit('third_coavailability_audit.json')
        self.check('third/descriptive_replication_status', a['n_libraries'] == 6 and
                   a['author_reported_biological_replicates_per_condition'] == 3 and
                   not a['cohort_p_values_computed'] and a['does_not_enter_original_200_label_null'])
        self.check('third/aggregation_audits', len(a['aggregation_checks']) == 6 and
                   all(x['matrix_gene_major_unique'] and x['all_assigned_counts_conserved'] for x in a['aggregation_checks']))
        m = self.read('third_library_metadata.tsv')
        self.check('third/library_unit', m['sample'].is_unique and m.condition.value_counts().to_dict() == {'Sham': 3, 'MCAO': 3})
        scores, effects = self.read('third_all_covered_coavailability_library_scores.tsv.gz'), self.read('third_all_covered_coavailability_effects.tsv.gz')
        self.check('third/unique_library_edge', not scores.duplicated(KEY + ['sample']).any())
        mx = scores.pivot(index=KEY, columns='sample', values='score').reindex(columns=m['sample'])
        case = m.condition.to_numpy() == 'MCAO'
        delta = mx.to_numpy()[:, case].mean(axis=1) - mx.to_numpy()[:, ~case].mean(axis=1)
        self.check('third/equal_library_effects', np.isfinite(mx).all().all() and
                   np.allclose(effects.set_index(KEY).loc[mx.index].score_difference, delta, atol=TOL, rtol=0))
        pairs = self.read('third_cross_cohort_concordance_edges.tsv.gz')
        s = self.read('third_cross_cohort_concordance_summary.tsv')
        for r in s.itertuples(index=False):
            d = pairs[pairs.baseline_dataset == r.baseline_dataset]
            av = d.score_difference_baseline.to_numpy(copy=True); bv = d.score_difference_third.to_numpy(copy=True)
            av[np.abs(av) <= TOL] = 0; bv[np.abs(bv) <= TOL] = 0
            nz = (av != 0) & (bv != 0); same = nz & (np.sign(av) == np.sign(bv))
            self.check('third/' + r.baseline_dataset + '/concordance', len(d) == r.n_common_eligible_edges and
                       same.sum() == r.n_same_direction and np.isclose(same.mean(), r.same_direction_fraction, atol=TOL) and
                       np.isclose(stats.spearmanr(av, bv).statistic, r.spearman_rho, atol=TOL))

    def optional(self):
        """Completion is audited model/output availability, not assumed convergence."""
        p = WORK / 'repro_decontx/completion_status.json'
        status = json.loads(p.read_text()) if p.exists() else {}
        needed = [TABLES / 'ambient_summary_audit.json', TABLES / 'rank_diagnostic_decontx_audit.json',
                  TABLES / 'ambient_rank_matched_audit.json']
        selection_pending = []
        if status.get('completed_samples') == 11:
            from repro_decontx_source import resolve
            for row in status.get('samples', []):
                try:
                    resolve(row['dataset'], row['sample'])
                except (AssertionError, FileNotFoundError) as exc:
                    selection_pending.append(str(exc))
        if status.get('completed_samples') != 11 or not all(x.exists() for x in needed) or selection_pending:
            self.pending.append({'stage': 'ambient', 'reason': 'Requires 11 convergence-rule-selected model audits, custom and matched summaries, and 22 corrected-network audits.',
                                 'pending_numerical_extensions': selection_pending,
                                 'completed_model_libraries': status.get('completed_samples', 0)})
        else:
            self.stage('ambient_completed_outputs', self.ambient)
        p = TABLES / 'third_rank_diagnostic_audit.json'
        if not p.exists():
            self.pending.append({'stage': 'third_liana', 'reason': 'Six completed third-cohort full-network audits and final descriptive summary are not yet all available.'})
        else:
            self.stage('third_liana_completed_outputs', self.third_liana)

    def ambient(self):
        a = json.loads((WORK / 'repro_decontx/completion_status.json').read_text())
        self.check('ambient/completed_status', a['expected_samples'] == a['completed_samples'] == 11 and len(a['samples']) == 11)
        self.check('ambient/matrix_audit_status', a['parameters_versions_and_matrix_audits_checked'] and a['n_zero_corrected_total'] == 0)
        from repro_decontx_source import resolve
        selected = self.read('ambient_selected_sources.tsv')
        self.check('ambient/11_selected_sources', len(selected) == 11 and not selected.duplicated(['dataset', 'sample']).any())
        nonconverged = []
        resolved = {}
        for row in a['samples']:
            p = WORK / 'repro_decontx' / row['dataset'] / row['sample'] / 'audit.json'
            rec = json.loads(p.read_text())
            self.check('ambient/' + row['sample'] + '/verified', all(rec['python_verification'].values()) and
                       rec['requested_parameters']['seed'] == 20260911 and rec['n_zero_corrected_total'] == 0)
            folder, decision = resolve(row['dataset'], row['sample'])
            chosen = json.loads((folder / 'audit.json').read_text())
            resolved[row['dataset'], row['sample']] = decision
            listed = selected[(selected.dataset == row['dataset']) & (selected['sample'] == row['sample'])]
            self.check('ambient/' + row['sample'] + '/convergence_source_hashes', len(listed) == 1 and
                       listed.selected_audit_sha256.iloc[0] == decision['selected_audit_sha256'] and
                       listed.default_audit_sha256.iloc[0] == decision['default_audit_sha256'] and
                       listed.selected_relative_directory.iloc[0] == decision['selected_relative_directory'])
            self.check('ambient/' + row['sample'] + '/selected_matrix_verified', all(chosen['python_verification'].values()))
            if not decision['selected_converged']:
                nonconverged.append(row['sample'])
        if nonconverged:
            self.warnings.append({'stage': 'ambient', 'nonconverged_samples': nonconverged,
                                  'meaning': 'Finished model runs are not equivalent to convergence; disclose these statuses and sensitivity decisions.'})
        a = self.audit('ambient_summary_audit.json')
        self.check('ambient/descriptive', a['n_libraries'] == 11 and a['no_new_pvalues'])
        for name in ['ambient_coverage.tsv', 'ambient_concordance_summary.tsv', 'ambient_fixed_examples.tsv',
                     'ambient_cross_cohort_effects.tsv.gz', 'ambient_markers_per_animal.tsv']:
            d = self.read(name)
            self.check(name + '/no_infinities', not np.isinf(d.select_dtypes(include=np.number)).any().any())
        self.diagnostic_networks('decontx')
        networks = json.loads((TABLES / 'rank_diagnostic_decontx_audit.json').read_text())
        self.check('ambient/network_source_selection_matches_resolver',
                   all(x['corrected_source_selection'] == resolved[x['dataset'], x['sample']] for x in networks['audits']))
        matched = self.audit('ambient_rank_matched_audit.json')
        self.check('ambient/matched_no_new_inference', not matched['new_pvalues'] and matched['raw_and_corrected_fixed_universe'] and
                   matched['n_config_library_combinations'] == 22 and len(matched['metrics']) == 6)
        self.check('ambient/matched_script_hash', matched['source_script_sha256'] == sha(WORK / 'repro_ambient_rank_comparison.py'))
        keys = self.read('ambient_rank_matched_candidate_keys.tsv')
        scores = self.read('ambient_rank_matched_animal_scores.tsv.gz')
        effects = self.read('ambient_rank_matched_effects.tsv.gz')
        summary = self.read('ambient_rank_matched_summary.tsv')
        for (config, source, acc, metric), block in scores.groupby(['config', 'count_source', 'dataset', 'metric']):
            m = block[['sample', 'condition']].drop_duplicates().sort_values('sample')
            self.check('ambient/matched_animal_unit', m['sample'].is_unique and m.condition.value_counts().to_dict() == GROUPS[acc])
            x = block.pivot(index=KEY, columns='sample', values='value').reindex(columns=m['sample'])
            wanted = pd.MultiIndex.from_frame(keys[keys.config == config][KEY])
            self.check('ambient/matched_common_finite', len(x) == len(wanted) and set(x.index) == set(wanted) and np.isfinite(x).all().all())
            case = m.condition.eq('MCAO').to_numpy()
            delta = x.to_numpy()[:, case].mean(axis=1) - x.to_numpy()[:, ~case].mean(axis=1)
            delta[np.abs(delta) <= TOL] = 0
            e = effects[(effects.config == config) & (effects.count_source == source) & (effects.dataset == acc) & (effects.metric == metric)].set_index(KEY).reindex(x.index)
            self.check('ambient/matched_equal_animal_effects', np.allclose(e.mcao_minus_sham, delta, atol=TOL, rtol=0))
        for row in summary.itertuples(index=False):
            self.check('ambient/matched_denominators', row.n_same_direction <= row.n_nonzero_both <= row.n_candidates and
                       row.n_candidates == len(keys[keys.config == row.config]) and
                       np.isclose(row.same_direction_all_fraction, row.n_same_direction / row.n_candidates, atol=TOL) and
                       (row.n_nonzero_both == 0 or np.isclose(row.same_direction_nonzero_fraction, row.n_same_direction / row.n_nonzero_both, atol=TOL)))

    def third_liana(self):
        a = self.audit('third_rank_diagnostic_audit.json', 'repro_third_liana.py')
        self.check('third_liana/descriptive', a['n_third_libraries'] == 6 and not a['inferential_pvalues_computed'] and a['does_not_enter_original_200_label_null'])
        for value, expected in a['input_target_sha256'].items():
            path = resolve_recorded_path(value)
            self.check('third_liana/target_hash', path.exists() and sha(path) == expected)
        for name in ['summary.tsv', 'effects.tsv.gz', 'library_scores.tsv.gz', 'ties.tsv']:
            d = self.read('third_rank_diagnostic_' + name)
            self.check(name + '/no_infinities', not np.isinf(d.select_dtypes(include=np.number)).any().any())
        folder = WORK / 'repro_third_cohort/results/liana'
        # Discover the archived location from the input target paths; never load a matrix.
        third_targets = [resolve_recorded_path(x) for x in a['input_target_sha256'] if 'repro_third_liana/' in x.replace('\\', '/')]
        self.check('third_liana/six_target_files', len(third_targets) == 6)
        for target in third_targets:
            stem = str(target).replace('__target.tsv.gz', '')
            rec = json.loads(Path(stem + '.json').read_text())
            self.check('third_liana/full_network_audit', rec['status'] == 'complete' and rec['computed_all_cell_type_pairs'] and
                       sha(Path(stem + '__full_network.tsv.gz')) == rec['full_network_sha256'])


def validate(require_complete=False, report_dir=None):
    v = Validator()
    for name, fn in [('required_outputs', v.inventory), ('animal_nulls', v.nulls), ('conditional_null', v.conditional),
                     ('rank_and_background', v.ranks), ('raw_network_provenance', v.diagnostic_networks),
                     ('sorted_pool_effects', v.sorted), ('third_custom', v.third)]:
        v.stage(name, fn)
    v.optional()
    failed = [x for x in v.checks if not x['passed']]
    code = 1 if failed else 2 if require_complete and v.pending else 0
    report = {'status': 'failed' if failed else 'incomplete' if require_complete and v.pending else 'completed_stages_pass',
              'require_complete': require_complete, 'exit_code': code, 'read_only_sources': True,
              'no_cell_matrices_loaded': True, 'no_new_model_fits': True,
              'checks': v.checks, 'stages': v.stages, 'n_checks': len(v.checks), 'n_failed': len(failed),
              'pending_stages': v.pending, 'warnings': v.warnings, 'checked_result_sha256': v.files,
              'script_sha256': sha(Path(__file__)),
              'versions': {x: importlib.metadata.version(x) for x in ['numpy', 'pandas', 'scipy']},
              'limits': ['Recorded biological provenance is checked; library independence is not inferred from file counts.',
                         'Statistical/source integrity does not certify biological mechanisms or replace visual figure review.',
                         'Large raw/corrected cell matrices are outside this cached validation boundary.']}
    if report_dir:
        report_dir = Path(report_dir)
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / 'repro_validation_audit.json').write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
        pd.DataFrame(v.checks).to_csv(report_dir / 'repro_validation_checks.tsv', sep='\t', index=False)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--require-complete', action='store_true')
    ap.add_argument('--report-dir', type=Path, default=BASE / 'outputs/reproducibility_v2/validation')
    args = ap.parse_args()
    r = validate(args.require_complete, args.report_dir)
    print(json.dumps({k: r[k] for k in ['status', 'exit_code', 'n_checks', 'n_failed', 'stages', 'pending_stages', 'warnings']}, indent=2))
    for check in r['checks']:
        if not check['passed']:
            print('FAIL', check['check'], check['detail'])
    return r['exit_code']


if __name__ == '__main__':
    sys.exit(main())
