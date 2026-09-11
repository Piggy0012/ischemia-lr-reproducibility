"""Descriptive global-universe RRA intervention using only saved small tables."""
from __future__ import annotations
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
KEY = ['source', 'target', 'ligand_complex', 'receptor_complex']
RENAMES = dict(zip(KEY, ['sender', 'receiver', 'ligand', 'receptor']))
INPUTS = {}
CHECKS = []
TOL = 1e-12

def register(p):
    INPUTS[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p

def check(name, condition):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)

def main():
    register(ROOT / 'repro_rank_background_plan.md')
    pkg = Path(importlib.util.find_spec('liana').origin).parent
    source = register(pkg / 'method/_pipe_utils/_aggregate.py')
    spec = importlib.util.spec_from_file_location('official_rank_background_audit', source)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    fixed_keys = pd.read_csv(register(OUT / 'null_fixed_candidate_keys.tsv'), sep='\t').rename(columns={v:k for k,v in RENAMES.items()})
    output_rows, effect_rows, summary_rows, global_keys, redundancy_rows, contexts, zero_rows = [], [], [], [], [], [], []
    for config in ['primary', 'reference_singlet']:
        target = pd.MultiIndex.from_frame(fixed_keys[fixed_keys.config == config][KEY]).sort_values()
        libraries = []
        common = None
        for acc in ['GSE174574', 'GSE245386']:
            folder = ROOT / 'repro_liana_rank_diagnostic/raw' / acc
            paths = sorted(folder.glob(config + '__*__full_network.tsv.gz'))
            check(config + '/' + acc + '/expected_libraries', len(paths) == (6 if acc == 'GSE174574' else 5))
            for path in paths:
                d = pd.read_csv(register(path), sep='\t', usecols=KEY + ['sample', 'condition', 'lr_means', 'expr_prod', 'lrscore',
                    'magnitude_rank', 'diagnostic_unique_column_magnitude_rank']).set_index(KEY)
                # Original LIANA magnitude columns were float32. Restore the
                # serialized columns to that precision before rank aggregation;
                # rank normalization otherwise changes at approximately 1e-7.
                # Every original/native and diagnostic RRA is validated below.
                d[['lr_means','expr_prod','lrscore']] = d[['lr_means','expr_prod','lrscore']].astype(np.float32)
                check(str(path) + '/unique_keys', not d.index.duplicated().any())
                finite = np.isfinite(d[['lr_means', 'expr_prod', 'lrscore']].to_numpy()).all(axis=1)
                available = d.index[finite]
                common = available if common is None else common.intersection(available)
                sample = d['sample'].iloc[0]; condition = d.condition.iloc[0]
                ap = folder / (config + '__' + sample + '.json')
                audit = json.loads(register(ap).read_text(encoding='utf-8'))
                libraries.append((acc, sample, condition, d, audit))
        common = common.sort_values()
        check(config + '/fixed_global_includes_previous_targets', set(target) <= set(common))
        g = common.to_frame(index=False); g.insert(0, 'config', config)
        global_keys.append(g.rename(columns=RENAMES))
        for acc, sample, condition, d, audit in libraries:
            native_specs = audit['native_magnitude_specs']; unique_specs = audit['diagnostic_unique_score_specs']
            for col, specs in [('magnitude_rank', native_specs), ('diagnostic_unique_column_magnitude_rank', unique_specs)]:
                reconstructed = mod._rank_aggregate(d.reset_index(), specs, 'rra')
                error = float(np.max(np.abs(reconstructed - d[col].to_numpy())))
                check(f'{config}/{acc}/{sample}/{col}/roundtrip_rra_reconstruction', error < TOL)
            fixed = d.loc[common].copy()
            rank = mod._rank_aggregate(fixed.reset_index(), unique_specs, 'rra')
            fixed['fixed_global_unique_magnitude_priority'] = 1 - rank
            for universe, frame in [('original_complete_full_network', d), ('fixed_global_network_intersection', fixed)]:
                a = stats.rankdata(frame.expr_prod, method='average')
                b = stats.rankdata(frame.lrscore, method='average')
                redundancy_rows.append({'config': config, 'dataset': acc, 'sample': sample, 'universe': universe,
                    'n_edges': len(frame), 'rank_vectors_exactly_equal': bool(np.array_equal(a, b)),
                    'n_different_average_tie_ranks': int((a != b).sum()),
                    'rank_spearman_rho': float(stats.spearmanr(a, b).statistic),
                    'max_abs_rank_difference': float(np.max(np.abs(a-b)))})
            contexts.append({'config': config, 'dataset': acc, 'sample': sample, 'condition': condition,
                'n_original_full_network_edges': len(d), 'n_fixed_global_edges': len(common),
                'n_fixed_target_edges': len(target), 'n_fixed_global_celltype_pairs': len(common.to_frame(index=False)[['source','target']].drop_duplicates()),
                'n_fixed_global_celltypes': len(set(common.get_level_values('source')) | set(common.get_level_values('target')))})
            for metric, values in [('native_full_network_priority', 1 - d.loc[target].magnitude_rank.to_numpy()),
                ('unique_full_network_priority', 1 - d.loc[target].diagnostic_unique_column_magnitude_rank.to_numpy()),
                ('unique_fixed_global_priority', fixed.loc[target].fixed_global_unique_magnitude_priority.to_numpy())]:
                r = target.to_frame(index=False).rename(columns=RENAMES)
                r.insert(0, 'metric', metric); r.insert(0, 'condition', condition); r.insert(0, 'sample', sample)
                r.insert(0, 'dataset', acc); r.insert(0, 'config', config); r['value'] = values
                output_rows.append(r)
    animal = pd.concat(output_rows, ignore_index=True)
    outkey = list(RENAMES.values())
    for (config, metric), block in animal.groupby(['config', 'metric'], sort=False):
        effects = []
        for acc in ['GSE174574', 'GSE245386']:
            d = block[block.dataset == acc]
            m = d[['sample', 'condition']].drop_duplicates().sort_values('sample')
            x = d.pivot(index=outkey, columns='sample', values='value').reindex(columns=m['sample'])
            case = m.condition.to_numpy() == 'MCAO'
            val = x.to_numpy(); eff = val[:,case].mean(axis=1)-val[:,~case].mean(axis=1)
            eff[np.abs(eff) <= TOL] = 0; effects.append(eff)
            zero_all = (np.abs(val) <= TOL).all(axis=1)
            zero_rows.append({'config':config,'dataset':acc,'metric':metric,'n_candidates':len(eff),
                'n_zero_disease_effect':int((eff == 0).sum()),
                'n_priority_zero_in_all_animals':int(zero_all.sum()),
                'n_zero_effect_with_all_animals_at_priority_zero':int(((eff == 0) & zero_all).sum())})
            row = x.index.to_frame(index=False); row.insert(0, 'metric', metric)
            row.insert(0, 'dataset', acc); row.insert(0, 'config', config); row['mcao_minus_sham'] = eff
            effect_rows.append(row)
        a,b = effects; nz = (a != 0) & (b != 0); same = nz & (np.sign(a) == np.sign(b))
        summary_rows.append({'config': config, 'metric': metric, 'n_candidates': len(a), 'n_nonzero_both': int(nz.sum()),
            'n_same_direction': int(same.sum()), 'n_zero_discovery': int((a==0).sum()), 'n_zero_validation': int((b==0).sum()),
            'same_direction_all_fraction': float(same.mean()), 'same_direction_nonzero_fraction': float(same.sum()/nz.sum()),
            'spearman_rho': float(stats.spearmanr(a,b).statistic)})
    tables = {'summary.tsv':pd.DataFrame(summary_rows), 'animal_scores.tsv.gz':animal,
              'effects.tsv.gz':pd.concat(effect_rows,ignore_index=True), 'global_keys.tsv.gz':pd.concat(global_keys,ignore_index=True),
              'monotonic_redundancy.tsv':pd.DataFrame(redundancy_rows), 'context.tsv':pd.DataFrame(contexts),
              'zero_effects.tsv':pd.DataFrame(zero_rows)}
    for suffix, frame in tables.items():
        frame.to_csv(OUT / ('rank_background_' + suffix), sep='\t', index=False)
    meta = {'status':'complete', 'n_checks':len(CHECKS), 'checks':CHECKS, 'input_sha256':INPUTS,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'new_pvalues':False,
        'upstream_celltype_expression_recomputed':False,
        'magnitude_input_dtype':'float32 restored from TSV; both saved full-network RRA outputs reproduced for every library',
        'intervention':'same unique-column RRA on a fixed whole-network edge intersection, then same target extraction'}
    (OUT/'rank_background_audit.json').write_text(json.dumps(meta,indent=2,allow_nan=False),encoding='utf-8')
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print(pd.DataFrame(contexts).groupby('config').agg(n_fixed_global_edges=('n_fixed_global_edges','first'),n_fixed_global_celltype_pairs=('n_fixed_global_celltype_pairs','first')).to_string())
    print(pd.DataFrame(redundancy_rows).groupby('universe').agg(all_equal=('rank_vectors_exactly_equal','all'),max_different=('n_different_average_tie_ranks','max')).to_string())

if __name__ == '__main__':
    main()
