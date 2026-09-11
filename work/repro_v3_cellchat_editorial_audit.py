"""Read-only independent checks of finalized CellChat descriptive tables."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / 'repro_v3_cellchat/tables'
KEY = ['source', 'target', 'ligand', 'receptor']
METRICS = ['custom_coavailability', 'cellchat_probability',
           'cellchat_strength_percentile', 'LIANA_unique_column_RRA_priority']


def read(name):
    return pd.read_csv(TABLES / name, sep='\t', float_precision='round_trip')


def file_record(path):
    path = Path(path)
    return {'path': path.relative_to(ROOT.parent).as_posix(),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    sources = ['summary_audit.json', 'target_sample_scores.tsv.gz',
               'target_disease_effects.tsv.gz', 'cross_cohort_concordance.tsv',
               'zero_strength_target_audit.tsv', 'tie_zero_audit.tsv',
               'fixed_complete_network_keys.tsv', 'fixed_complete_target_keys.tsv',
               'sample_manifest.tsv', 'original_187_to_controlled_32_attrition.json']
    audit = json.loads((TABLES / sources[0]).read_text())
    assert audit['status'] == 'complete'
    scores = read('target_sample_scores.tsv.gz')
    target = read('fixed_complete_target_keys.tsv')
    fixed = read('fixed_complete_network_keys.tsv')
    manifest = read('sample_manifest.tsv')
    summary = read('cross_cohort_concordance.tsv')
    effects = read('target_disease_effects.tsv.gz')
    zeros = read('zero_strength_target_audit.tsv').set_index(KEY)
    ties = read('tie_zero_audit.tsv').set_index('sample')
    target_index = pd.MultiIndex.from_frame(target)
    assert len(fixed) == 1222 and len(target) == 32 and len(manifest) == 11
    assert len(target[['ligand', 'receptor']].drop_duplicates()) == 29
    assert target.groupby(['source', 'target']).size().to_dict() == {
        ('Astrocyte', 'Endothelial'): 8, ('Endothelial', 'Astrocyte'): 24}
    assert not scores.duplicated(KEY + ['sample']).any()
    assert len(scores) == 352
    checked = []
    epsilon_changes = []
    for metric in METRICS:
        vectors = []
        for acc in ['GSE174574', 'GSE245386']:
            d = scores[scores.dataset == acc]
            # Aggregate directly by condition rather than reusing the saved
            # contrast calculation or relying on manifest column order.
            means = d.groupby(KEY + ['condition'], sort=True)[metric].mean().unstack('condition')
            raw = (means.MCAO - means.Sham).reindex(target_index).to_numpy()
            contrast = np.where(np.abs(raw) <= 1e-12, 0.0, raw)
            saved = effects[(effects.dataset == acc) & (effects.metric == metric)].set_index(KEY).loc[target_index]
            assert np.allclose(raw, saved.raw_mcao_minus_sham, rtol=1e-12, atol=1e-16)
            assert np.allclose(contrast, saved.mcao_minus_sham, rtol=1e-12, atol=1e-16)
            vectors.append(contrast)
            epsilon_changes.append({'metric': metric, 'dataset': acc,
                'nonzero_contrasts_zeroed_at_1e_12': int(((raw != 0) & (contrast == 0)).sum())})
        left, right = vectors
        nonzero = (left != 0) & (right != 0)
        same = nonzero & (np.sign(left) == np.sign(right))
        rho = float(spearmanr(left, right).statistic)
        saved = summary[(summary.scope == 'all_fixed_targets') & (summary.metric == metric)].iloc[0]
        assert saved.n_candidates == 32
        assert saved.n_nonzero_both == nonzero.sum()
        assert saved.n_same_direction == same.sum()
        assert saved.n_zero_discovery == (left == 0).sum()
        assert saved.n_zero_external == (right == 0).sum()
        assert np.isclose(rho, saved.spearman_rho, rtol=1e-12, atol=1e-14)
        checked.append({'metric': metric, 'n_candidates': 32,
            'n_nonzero_both': int(nonzero.sum()), 'n_same_direction': int(same.sum()),
            'same_direction_all_fraction': float(same.sum()/32),
            'same_direction_nonzero_fraction': float(same.sum()/nonzero.sum()) if nonzero.any() else None,
            'n_zero_discovery': int((left == 0).sum()), 'n_zero_external': int((right == 0).sum()),
            'spearman_rho': rho})
    zero_flags = scores.groupby(KEY).cellchat_probability.apply(lambda x: bool((x == 0).all())).reindex(target_index)
    assert np.array_equal(zero_flags.to_numpy(), zeros.loc[target_index].all_11_strength_zero.to_numpy())
    zero_targets = target_index[zero_flags.to_numpy()]
    rank_e = effects[effects.metric == 'cellchat_strength_percentile'].set_index(KEY)
    allzero_rank_changes = int(rank_e.loc[zero_targets].groupby(level=KEY).mcao_minus_sham.apply(lambda x: bool((x != 0).any())).sum())
    assert allzero_rank_changes == audit['all_11_zero_strength_targets_with_nonzero_derived_disease_effect']
    assert int(zero_flags.sum()) == audit['all_11_zero_strength_targets']
    for sample, row in ties.iterrows():
        full = read(sample + '__fixed_full_network_scores.tsv.gz')
        strength = full.cellchat_probability.to_numpy()
        count = int((strength == 0).sum())
        assert count == row.cellchat_zero_network_rows
        # Closed-form average-rank percentile for the zero tie block.
        floor = (count - 1) / (2 * (len(fixed) - 1)) if count else None
        if count:
            assert np.allclose(full.loc[strength == 0, 'cellchat_strength_percentile'], floor, rtol=0, atol=1e-16)
            assert np.isclose(row.cellchat_zero_strength_percentile_floor, floor, rtol=0, atol=1e-16)
    result = {'status': 'complete', 'reviewed_at_utc': datetime.now(timezone.utc).isoformat(),
        'checks': ['fixed universe and direction counts', 'condition contrasts independently regrouped',
                   'same-direction counts and both denominators', 'Spearman correlations',
                   'all-library zero-strength target count', 'changing ranks of all-zero strengths',
                   'closed-form zero tie-block percentile floor'],
        'primary_metrics': checked, 'epsilon_changes': epsilon_changes,
        'all_11_zero_strength_targets': int(zero_flags.sum()),
        'all_11_zero_strength_targets_with_nonzero_rank_effect': allzero_rank_changes,
        'native_model_rerun': False, 'upstream_files_modified': False,
        'source_files': [file_record(TABLES / x) for x in sources],
        'review_script': file_record(__file__)}
    dest = ROOT / 'repro_v3_cellchat_editorial_audit.json'
    dest.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
