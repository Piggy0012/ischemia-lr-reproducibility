"""Gene-major, context-only cache and descriptive LIANA audit for GSE332910.

Run --stage cache first; --stage liana is intentionally a separate explicit step.
The six deposited libraries have unverified individual-animal provenance.
"""
from pathlib import Path
import argparse, gc, hashlib, importlib.metadata, json, os, time
for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ[name] = '1'
import numpy as np
import pandas as pd
import psutil
from scipy import sparse, stats
from stream_annotate import chunks

ROOT = Path(__file__).resolve().parent
COHORT = ROOT / 'repro_third_cohort'
PROCESSED = COHORT / 'processed/GSE332910'
CACHE = COHORT / 'context_cache'
OUT = ROOT / 'repro_third_liana'
TABLES = ROOT.parent / 'outputs/reproducibility_v2/tables'
ACC = 'GSE332910'
NATIVE_KEY = ['source', 'target', 'ligand_complex', 'receptor_complex']
KEY = ['sender', 'receiver', 'ligand', 'receptor']
RENAMES = dict(zip(NATIVE_KEY, KEY))
METRICS = ['custom_score', 'lr_means', 'expr_prod', 'lrscore',
           'native_magnitude_priority', 'diagnostic_unique_column_priority']
TOL = 1e-12

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def memory():
    info = psutil.Process().memory_info()
    return {'rss_bytes': info.rss, 'peak_working_set_bytes': getattr(info, 'peak_wset', info.rss),
            'system_available_bytes': psutil.virtual_memory().available}

def manifest():
    records = sorted(json.loads((COHORT / 'download_manifest.json').read_text(encoding='utf-8')), key=lambda x: x['sample'])
    assert len(records) == 6
    assert pd.Series([r['condition'] for r in records]).value_counts().to_dict() == {'MCAO': 3, 'Sham': 3}
    return records

def build_cache(rec):
    sample = rec['sample']
    folder = CACHE / sample
    folder.mkdir(parents=True, exist_ok=True)
    ap = folder / 'audit.json'
    if ap.exists():
        report = json.loads(ap.read_text())
        assert report['status'] == 'complete' and report['sample'] == sample
        print('CACHE EXISTS', sample, report['matrix_shape'], flush=True)
        return report
    started = time.monotonic()
    qp = PROCESSED / (sample + '_cells.tsv.gz')
    gp = PROCESSED / (sample + '_features.tsv.gz')
    q = pd.read_csv(qp, sep='\t', float_precision='round_trip')
    genes = pd.read_csv(gp, sep='\t', float_precision='round_trip')
    assigned = q.cell_type != 'Unassigned'
    counts = q.loc[assigned, 'cell_type'].value_counts()
    retained = counts[counts >= 30].index
    keep = assigned & q.cell_type.isin(retained)
    assert q.loc[keep, 'qc_pass'].all()
    assert set(['Astrocyte', 'Endothelial']) <= set(retained)
    sub = q.loc[keep].copy()
    sub.insert(0, 'raw_matrix_cell_index', np.flatnonzero(keep))
    sub['qualified_barcode'] = sample + ':' + sub.barcode.astype(str)
    assert sub.qualified_barcode.is_unique and genes.symbol.is_unique
    lookup = np.full(len(q), -1, np.int32)
    lookup[keep] = np.arange(len(sub), dtype=np.int32)
    expected_nnz = int(sub.n_genes.sum())
    ng = len(genes)
    # Source coordinates are gene-major, exactly CSC column order for cells x genes.
    # Only retained context entries are allocated; no all-nucleus COO/CSR is built.
    data = np.empty(expected_nnz, np.int32)
    indices = np.empty(expected_nnz, np.int32)
    per_gene = np.zeros(ng, np.int64)
    matrix = COHORT / next(x['destination_relative'] for x in rec['extracted_files'] if x['source_member'].endswith('matrix.mtx.gz'))
    it = chunks(matrix)
    dims = next(it)
    assert dims == rec['matrix_shape_genes_cells_nnz']
    assert dims[:2] == [ng, len(q)]
    seen = offset = 0
    previous = (-1, -1)
    tracked = [memory()]
    for block in it:
        gi, ci, x = block.T
        assert np.all(x > 0)
        assert (int(gi[0]), int(ci[0])) > previous
        assert np.all((np.diff(gi) > 0) | ((np.diff(gi) == 0) & (np.diff(ci) > 0)))
        previous = (int(gi[-1]), int(ci[-1]))
        seen += len(block)
        target_ci = lookup[ci]
        valid = target_ci >= 0
        k = int(valid.sum())
        data[offset:offset+k] = x[valid]
        indices[offset:offset+k] = target_ci[valid]
        per_gene += np.bincount(gi[valid], minlength=ng)
        offset += k
    assert seen == dims[2] and offset == expected_nnz
    indptr = np.r_[0, per_gene.cumsum()].astype(np.int32)
    x = sparse.csc_matrix((data, indices, indptr), shape=(len(sub), ng), copy=False)
    assert x.has_canonical_format
    tracked.append(memory())
    y = x.tocsr()
    del x, data, indices, indptr, per_gene, block, q
    gc.collect()
    assert y.has_canonical_format
    assert np.array_equal(np.asarray(y.sum(axis=1)).ravel(), sub.n_umis.to_numpy())
    assert np.array_equal(np.diff(y.indptr), sub.n_genes.to_numpy())
    tracked.append(memory())
    cp, sp, vp = folder / 'counts.npz', folder / 'cells.tsv.gz', folder / 'genes.tsv.gz'
    sparse.save_npz(cp, y, compressed=True)
    sub.to_csv(sp, sep='\t', index=False)
    genes.to_csv(vp, sep='\t', index=False)
    tracked.append(memory())
    report = {'status': 'complete', 'dataset': ACC, 'sample': sample, 'condition': rec['condition'],
        'selection': 'original primary labels; QC-passing assigned classes with >=30 cells per library',
        'matrix_shape': list(y.shape), 'nnz': y.nnz,
        'csr_array_bytes': int(y.data.nbytes + y.indices.nbytes + y.indptr.nbytes),
        'matrix_format': 'CSR cells x submitted gene symbols; raw integer UMI',
        'context_counts': sub.cell_type.value_counts().to_dict(),
        'source_raw_shape': dims, 'all_raw_entries_streamed': seen,
        'source_gene_major_unique': True, 'total_umi_and_detected_genes_match_qc': True,
        'n_verified_independent_animals': None, 'input_qc_sha256': sha(qp), 'input_gene_sha256': sha(gp),
        'input_matrix_sha256': sha(matrix),
        'output_sha256': {p.name: sha(p) for p in [cp, sp, vp]},
        'memory_observations': tracked, 'script_sha256': sha(__file__),
        'elapsed_seconds': time.monotonic() - started}
    ap.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('CACHE DONE', sample, y.shape, 'nnz', y.nnz, 'CSR MB', round(report['csr_array_bytes']/1e6, 1), 'memory', tracked[-1], flush=True)
    del y, sub, genes
    gc.collect()
    return report

def load_context(sample):
    folder = CACHE / sample
    audit = json.loads((folder / 'audit.json').read_text())
    for name, expected in audit['output_sha256'].items():
        assert sha(folder / name) == expected
    q = pd.read_csv(folder / 'cells.tsv.gz', sep='\t', float_precision='round_trip')
    genes = pd.read_csv(folder / 'genes.tsv.gz', sep='\t', float_precision='round_trip').symbol.astype(str).to_numpy()
    x = sparse.load_npz(folder / 'counts.npz').tocsr()
    assert list(x.shape) == audit['matrix_shape']
    assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(), q.n_umis.to_numpy())
    prior = pd.read_csv(PROCESSED/(sample+'_pseudobulk.tsv.gz'), sep='\t', index_col=0, float_precision='round_trip')
    assert np.array_equal(prior.index.astype(str), genes)
    for typ in q.cell_type.unique():
        sums = np.asarray(x[(q.cell_type == typ).to_numpy()].sum(axis=0)).ravel()
        assert np.array_equal(sums, prior['primary__'+typ].to_numpy()), 'Context cache differs from independent full-stream pseudobulk'
    return x, q, genes, audit

def process(rec):
    # Import the full inference stack only at the separately authorized run stage.
    import anndata as ad
    import liana as li
    from liana.method._pipe_utils._aggregate import _rank_aggregate
    import liana.method._pipe_utils._aggregate as aggregate_module
    sample = rec['sample']
    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / ('primary__' + sample)
    ap = Path(str(stem) + '.json')
    fp = Path(str(stem) + '__full_network.tsv.gz')
    tp = Path(str(stem) + '__target.tsv.gz')
    if all(p.exists() for p in [ap, fp, tp]):
        report = json.loads(ap.read_text())
        assert report['status'] == 'complete'
        print('LIANA EXISTS', sample, flush=True)
        return report
    started = time.monotonic()
    x, q, genes, ca = load_context(sample)
    totals = np.asarray(x.sum(axis=1)).ravel().astype(float)
    y = x.astype(np.float32)
    del x
    gc.collect()
    y = y.multiply((1e4 / totals).astype(np.float32)[:, None]).tocsr()
    y.data = np.log1p(y.data)
    n_explicit_zeros = int((y.data == 0).sum())
    assert n_explicit_zeros == 0 and np.isfinite(y.data).all() and np.all(y.data > 0), 'Raw integer input normalization must retain strictly positive finite stored values'
    a = ad.AnnData(y, obs=pd.DataFrame({'cell_type': pd.Categorical(q.cell_type.to_numpy())},
                   index=q.qualified_barcode.astype(str)), var=pd.DataFrame(index=genes))
    rp = ROOT / 'literature/mouseconsensus.csv'
    resource = pd.read_csv(rp)[['source_genesymbol', 'target_genesymbol']].drop_duplicates().rename(
        columns={'source_genesymbol': 'ligand', 'target_genesymbol': 'receptor'})
    print('LIANA START', sample, a.shape, memory(), flush=True)
    result = li.mt.rank_aggregate(a, groupby='cell_type', resource=resource, expr_prop=.1, min_cells=30,
        use_raw=False, n_perms=100, seed=20260910, n_jobs=1, return_all_lrs=False, inplace=False, verbose=False)
    assert not result.duplicated(NATIVE_KEY).any()
    specs = li.mt.rank_aggregate.magnitude_specs
    unique, used = {}, set()
    for name, spec in specs.items():
        if spec[0] not in used:
            unique[name] = spec
            used.add(spec[0])
    native_error = float(np.max(np.abs(_rank_aggregate(result.copy(), specs, 'rra') - result.magnitude_rank.to_numpy())))
    assert native_error < TOL
    result['diagnostic_unique_column_magnitude_rank'] = _rank_aggregate(result.copy(), unique, 'rra')
    result['native_magnitude_priority'] = 1 - result.magnitude_rank
    result['diagnostic_unique_column_priority'] = 1 - result.diagnostic_unique_column_magnitude_rank
    for key, value in [('count_source', 'raw'), ('config', 'primary'), ('dataset', ACC),
                       ('sample', sample), ('condition', rec['condition']), ('expression_eligible', True)]:
        result[key] = value
    target = result.loc[((result.source == 'Astrocyte') & (result.target == 'Endothelial')) |
                        ((result.source == 'Endothelial') & (result.target == 'Astrocyte'))].copy()
    result.to_csv(fp, sep='\t', index=False)
    target.to_csv(tp, sep='\t', index=False)
    report = {'status': 'complete', 'dataset': ACC, 'sample': sample, 'condition': rec['condition'],
        'config': 'primary', 'count_source': 'raw', 'n_verified_independent_animals': None,
        'author_reported_biological_replicates_per_condition': 3,
        'liana_version': importlib.metadata.version('liana'),
        'methods': [m.method_name for m in li.mt.rank_aggregate.methods],
        'native_magnitude_specs': specs, 'diagnostic_unique_score_specs': unique,
        'parameters': {'expr_prop': .1, 'min_cells': 30, 'use_raw': False, 'n_perms': 100,
                       'seed': 20260910, 'n_jobs': 1, 'return_all_lrs': False, 'aggregate_method': 'rra'},
        'all_pair_rows': len(result), 'exported_target_rows': len(target),
        'computed_all_cell_type_pairs': True, 'context_counts': ca['context_counts'],
        'n_context_cells': len(q), 'n_context_genes': len(genes),
        'denominator': 'per-cell total raw UMI across every submitted gene, identical to original QC total',
        'normalization': 'float32 log1p(CP10k), same operations as original two-cohort rank diagnostic',
        'normalized_stored_values_strictly_positive': True,
        'normalized_explicit_zeros': n_explicit_zeros,
        'all_context_pseudobulk_matches_independent_stream_aggregation': True,
        'native_full_network_reconstruction_max_abs_error': native_error,
        'new_cohort_has_no_prior_liana_export': True, 'raw_target_scores_verified_against_prior': None,
        'package_modified': False, 'diagnostic_is_official_package_release': False,
        'resource_sha256': sha(rp), 'context_cache_audit_sha256': sha(CACHE/sample/'audit.json'),
        'full_network_sha256': sha(fp), 'target_sha256': sha(tp),
        'source_aggregate_sha256': sha(aggregate_module.__file__), 'script_sha256': sha(__file__),
        'plan_sha256': sha(ROOT/'repro_third_liana_plan.md'),
        'memory_after_inference': memory(), 'elapsed_seconds': time.monotonic()-started}
    ap.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print('LIANA DONE', sample, len(result), len(target), report['memory_after_inference'], flush=True)
    del result, target, a, y, q, genes
    gc.collect()
    return report

def custom_scores(acc, m, common):
    folder = PROCESSED if acc == ACC else ROOT/'processed'/acc
    result = np.empty((len(common), len(m)), float)
    for j, sample in enumerate(m['sample']):
        mean = pd.read_csv(folder/(sample+'_mean_logcp10k.tsv.gz'), sep='\t', index_col=0, float_precision='round_trip')
        present = set(mean.index)
        for i, (sender, receiver, ligand, receptor) in enumerate(common):
            lg, rg = ligand.split('_'), receptor.split('_')
            assert set(lg+rg) <= present
            result[i, j] = np.sqrt(mean.loc[lg, 'primary__'+sender].min() * mean.loc[rg, 'primary__'+receiver].min())
    return pd.DataFrame(result, index=common, columns=m['sample'])

def summarize():
    datasets = ['GSE174574', 'GSE245386', ACC]
    targets, matrices, manifests, complete = {}, {}, {}, {}
    hashes = {}
    for acc in datasets:
        folder = OUT if acc == ACC else ROOT/'repro_liana_rank_diagnostic/raw'/acc
        paths = sorted(folder.glob('primary__*__target.tsv.gz'))
        assert len(paths) == (5 if acc == 'GSE245386' else 6), (acc, len(paths))
        targets[acc] = pd.concat([pd.read_csv(p, sep='\t', float_precision='round_trip').rename(columns=RENAMES) for p in paths], ignore_index=True)
        d = targets[acc]
        assert not d.duplicated(KEY+['sample']).any()
        assert set(d.dataset) == {acc} and set(d.config) == {'primary'}
        manifests[acc] = d[['sample', 'condition']].drop_duplicates().sort_values('sample').reset_index(drop=True)
        assert manifests[acc].condition.value_counts().to_dict() == {'Sham': 3, 'MCAO': 2 if acc == 'GSE245386' else 3}
        matrices[acc] = {}
        ok = None
        for metric in METRICS[1:]:
            mat = d.pivot(index=KEY, columns='sample', values=metric).reindex(columns=manifests[acc]['sample'])
            rows = mat.index[np.isfinite(mat.to_numpy()).all(axis=1)]
            ok = rows if ok is None else ok.intersection(rows)
            matrices[acc][metric] = mat
        complete[acc] = ok.sort_values()
        for p in paths:
            hashes[str(p.relative_to(ROOT))] = sha(p)
    global_common = complete[datasets[0]].intersection(complete[datasets[1]]).intersection(complete[ACC]).sort_values()
    comparisons = []
    for baseline in datasets[:2]:
        pair = complete[baseline].intersection(complete[ACC]).sort_values()
        comparisons.extend([(baseline, 'pairwise_complete', pair), (baseline, 'all_three_cohorts_complete', global_common)])
    summaries, effects, scores, ties = [], [], [], []
    for baseline, universe, common in comparisons:
        if not len(common):
            for metric in METRICS:
                summaries.append({'baseline_dataset':baseline, 'third_dataset':ACC, 'universe':universe,
                    'config':'primary', 'metric':metric, 'n_candidates':0, 'n_nonzero_both':0,
                    'n_same_direction':0, 'same_direction_all_fraction':None,
                    'same_direction_nonzero_fraction':None, 'n_zero_baseline':0, 'n_zero_third':0,
                    'spearman_rho':None, 'status':'no_complete_candidates', 'inferential_pvalues_computed':False})
            continue
        for acc in [baseline, ACC]:
            matrices[acc]['custom_score'] = custom_scores(acc, manifests[acc], common)
        for metric in METRICS:
            deltas = []
            for acc in [baseline, ACC]:
                mat = matrices[acc][metric].loc[common]
                case = manifests[acc].condition.to_numpy() == 'MCAO'
                x = mat.to_numpy()
                delta = x[:,case].mean(axis=1) - x[:,~case].mean(axis=1)
                delta[np.abs(delta) <= TOL] = 0
                deltas.append(delta)
                e = common.to_frame(index=False)
                e['mcao_minus_sham'] = delta
                e['dataset'], e['baseline_dataset'], e['universe'], e['metric'] = acc, baseline, universe, metric
                effects.append(e)
                s = mat.rename_axis(columns='sample').stack().rename('value').reset_index()
                s['dataset'], s['baseline_dataset'], s['universe'], s['metric'] = acc, baseline, universe, metric
                scores.append(s.merge(manifests[acc], on='sample', validate='many_to_one'))
            a, b = deltas
            nz = (a != 0) & (b != 0)
            same = nz & (np.sign(a) == np.sign(b))
            summaries.append({'baseline_dataset': baseline, 'third_dataset': ACC, 'universe': universe,
                'config': 'primary', 'metric': metric, 'n_candidates': len(common),
                'n_nonzero_both': int(nz.sum()), 'n_same_direction': int(same.sum()),
                'same_direction_all_fraction': float(same.mean()),
                'same_direction_nonzero_fraction': float(same.sum()/nz.sum()) if nz.sum() else None,
                'n_zero_baseline': int((a == 0).sum()), 'n_zero_third': int((b == 0).sum()),
                'spearman_rho': float(stats.spearmanr(a,b).statistic) if np.ptp(a) > 0 and np.ptp(b) > 0 else None,
                'status': 'complete' if np.ptp(a) > 0 and np.ptp(b) > 0 else 'rho_undefined_constant_effect',
                'inferential_pvalues_computed': False})
        for acc in [baseline, ACC]:
            for sample, d in targets[acc].groupby('sample'):
                selected = d.set_index(KEY).loc[common]
                for metric in ['magnitude_rank', 'diagnostic_unique_column_magnitude_rank']:
                    r = selected[metric].to_numpy()
                    c = np.unique(r, return_counts=True)[1]
                    ties.append({'baseline_dataset': baseline, 'dataset': acc, 'sample': sample,
                        'universe': universe, 'metric': metric, 'n_edges': len(r),
                        'n_rank_one': int((np.abs(r-1) <= TOL).sum()),
                        'rank_one_fraction': float((np.abs(r-1) <= TOL).mean()),
                        'n_edges_in_exact_tie_groups': int(c[c>1].sum()), 'largest_exact_tie_group': int(c.max())})
    TABLES.mkdir(parents=True, exist_ok=True)
    for name, frame in [('summary.tsv', pd.DataFrame(summaries)),
                        ('effects.tsv.gz', pd.concat(effects, ignore_index=True) if effects else pd.DataFrame(columns=KEY+['dataset','metric','mcao_minus_sham'])),
                        ('library_scores.tsv.gz', pd.concat(scores, ignore_index=True) if scores else pd.DataFrame(columns=KEY+['dataset','sample','metric','value'])),
                        ('ties.tsv', pd.DataFrame(ties))]:
        frame.to_csv(TABLES/('third_rank_diagnostic_'+name), sep='\t', index=False)
    audit = {'status': 'complete', 'third_dataset': ACC, 'n_third_libraries': 6,
        'n_verified_independent_third_animals': None, 'count_source': 'raw', 'config': 'primary',
        'author_reported_biological_replicates_per_condition': 3,
        'candidate_rule': 'All selected libraries and all LIANA metrics finite; identical edge universe for all six metric families within each comparison.',
        'pairwise_complete_library_counts': {'GSE174574_vs_GSE332910': 12, 'GSE245386_vs_GSE332910': 11},
        'all_three_cohorts_complete_library_count': 17, 'all_three_cohorts_complete_candidates': len(global_common),
        'no_missing_LIANA_scores_imputed': True, 'inferential_pvalues_computed': False,
        'tsv_float_parser': 'pandas float_precision=round_trip for all saved numerical inputs; preserves persisted double-precision ties',
        'does_not_enter_original_200_label_null': True, 'input_target_sha256': hashes,
        'script_sha256': sha(__file__), 'summary': summaries}
    (TABLES/'third_rank_diagnostic_audit.json').write_text(json.dumps(audit, indent=2, allow_nan=False), encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['cache', 'liana', 'summary'], default='cache')
    parser.add_argument('--samples', nargs='+')
    args = parser.parse_args()
    records = [r for r in manifest() if not args.samples or r['sample'] in args.samples]
    assert records
    if args.stage == 'summary':
        summarize()
    else:
        for rec in records:
            (build_cache if args.stage == 'cache' else process)(rec)
            if memory()['peak_working_set_bytes'] > 2 * 1024**3:
                raise RuntimeError('Observed process peak exceeds 2 GiB; stop before another library and coordinate available memory.')
