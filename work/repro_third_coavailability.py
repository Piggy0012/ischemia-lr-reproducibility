"""Streaming descriptive coavailability audit for GSE332910.

No animal-level P values: independence of the six submitted libraries is not
established. All resource-covered edges are exported before detection gating.
No absent feature rows are silently imputed for the third-cohort comparison.
"""
from pathlib import Path
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, hashlib, json
import numpy as np
import pandas as pd
from scipy import stats
from stream_annotate import MARKERS, chunks

BASE = Path(__file__).resolve().parent
ROOT = BASE / 'repro_third_cohort'
ACC = 'GSE332910'
FOLDER = ROOT / 'processed' / ACC
OUT = BASE.parent / 'outputs' / 'reproducibility_v2' / 'tables'
KEY = ['sender', 'receiver', 'ligand', 'receptor']
EXAMPLES = [
    ('Astrocyte', 'Endothelial', 'Spp1', 'Itga5_Itgb1'),
    ('Astrocyte', 'Endothelial', 'Timp3', 'Kdr'),
    ('Endothelial', 'Astrocyte', 'Ptn', 'Ptprz1'),
    ('Endothelial', 'Astrocyte', 'Plat', 'Lrp1'),
    ('Astrocyte', 'Endothelial', 'Col4a1', 'Itga3_Itgb1'),
]

def aggregate(rec):
    sample, condition = rec['sample'], rec['condition']
    marker = FOLDER / (sample + '_stream_aggregation.json')
    if marker.exists():
        return json.loads(marker.read_text())
    q = pd.read_csv(FOLDER / (sample + '_cells.tsv.gz'), sep='\t')
    g = pd.read_csv(FOLDER / (sample + '_features.tsv.gz'), sep='\t')
    assert g.symbol.is_unique and q.barcode.is_unique
    typ = list(MARKERS)
    tid = {x: i for i, x in enumerate(typ)}
    labels = q.cell_type.map(tid).fillna(-1).to_numpy(int)
    n = np.bincount(labels[labels >= 0], minlength=len(typ))
    total = q.n_umis.to_numpy(float)
    ng = len(g)
    umi = np.zeros((len(typ), ng), np.float64)
    logs = np.zeros_like(umi)
    detected = np.zeros((len(typ), ng), np.int64)
    matrix = ROOT / next(x['destination_relative'] for x in rec['extracted_files'] if x['source_member'].endswith('matrix.mtx.gz'))
    iterator = chunks(matrix)
    dims = next(iterator)
    assert dims == rec['matrix_shape_genes_cells_nnz']
    seen = 0
    previous = (-1, -1)
    for block in iterator:
        gi, ci, counts = block.T
        # Strict gene-major ordering also verifies that a cell-gene pair occurs once.
        assert np.all(counts > 0)
        assert (int(gi[0]), int(ci[0])) > previous
        assert np.all((np.diff(gi) > 0) | ((np.diff(gi) == 0) & (np.diff(ci) > 0)))
        previous = (int(gi[-1]), int(ci[-1]))
        seen += len(block)
        valid = labels[ci] >= 0
        gi, ci, counts = gi[valid], ci[valid], counts[valid]
        index = labels[ci] * ng + gi
        umi += np.bincount(index, weights=counts, minlength=ng * len(typ)).reshape(len(typ), ng)
        logs += np.bincount(index, weights=np.log1p(counts * 1e4 / total[ci]), minlength=ng * len(typ)).reshape(len(typ), ng)
        detected += np.bincount(index, minlength=ng * len(typ)).reshape(len(typ), ng)
    assert seen == dims[2]
    expected_umi = np.bincount(labels[labels >= 0], weights=total[labels >= 0], minlength=len(typ))
    assert np.array_equal(umi.sum(axis=1), expected_umi)
    assert np.all(detected <= n[:, None])
    mean = logs / np.maximum(n[:, None], 1)
    fraction = detected / np.maximum(n[:, None], 1)
    columns = ['primary__' + x for x in typ]
    for suffix, values in [('pseudobulk', umi.astype(np.int64)), ('mean_logcp10k', mean), ('fractions', fraction)]:
        d = pd.DataFrame(values.T, index=g.symbol.to_numpy(), columns=columns)
        d.index.name = 'gene'
        d.to_csv(FOLDER / (sample + '_' + suffix + '.tsv.gz'), sep='\t')
    rows = [{'dataset': ACC, 'sample': sample, 'condition': condition, 'config': 'primary', 'cell_type': t,
             'n_cells': int(n[j]), 'total_target_counts': int(umi[j].sum())} for j, t in enumerate(typ)]
    pd.DataFrame(rows).to_csv(FOLDER / (sample + '_sample_counts.tsv'), sep='\t', index=False)
    audit = {'sample': sample, 'processed_nonzero_entries': seen, 'matrix_gene_major_unique': True,
             'all_assigned_counts_conserved': True, 'n_target_astrocytes': int(n[tid['Astrocyte']]),
             'n_target_endothelial': int(n[tid['Endothelial']]), 'all_class_summaries': True}
    marker.write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print('AGGREGATED', sample, audit['n_target_astrocytes'], audit['n_target_endothelial'], flush=True)
    return audit

def compare(edges):
    summaries, comparisons = [], []
    for baseline in ['GSE174574', 'GSE245386']:
        prior = pd.read_csv(BASE / 'results' / baseline / 'primary__communication.tsv.gz', sep='\t')
        q = prior[prior.eligible_10][KEY + ['score_difference']].merge(edges[edges.eligible_10][KEY + ['score_difference']], on=KEY, suffixes=('_baseline', '_third'), validate='1:1')
        a, b = q.score_difference_baseline.to_numpy(copy=True), q.score_difference_third.to_numpy(copy=True)
        a[np.abs(a) <= 1e-12] = 0
        b[np.abs(b) <= 1e-12] = 0
        same = ((np.sign(a) == np.sign(b)) & (a != 0) & (b != 0))
        q['same_nonzero_direction'] = same
        q['baseline_dataset'] = baseline
        comparisons.append(q)
        summaries.append({'baseline_dataset': baseline, 'third_dataset': ACC, 'baseline_config': 'primary',
                          'third_config': 'primary', 'gate': '10% constituent detection in >=2 libraries of either condition, separately for ligand and receptor, in each cohort',
                          'n_common_eligible_edges': len(q), 'n_same_direction': int(same.sum()),
                          'same_direction_fraction': float(same.mean()), 'n_nonzero_both': int(((a != 0) & (b != 0)).sum()),
                          'spearman_rho': float(stats.spearmanr(a, b).statistic),
                          'inference': 'author-reported n=3 biological replicates/group; mice per library unspecified; descriptive, no P value'})
    pd.DataFrame(summaries).to_csv(OUT / 'third_cross_cohort_concordance_summary.tsv', sep='\t', index=False)
    pd.concat(comparisons, ignore_index=True).to_csv(OUT / 'third_cross_cohort_concordance_edges.tsv.gz', sep='\t', index=False)
    return summaries


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    records = sorted(json.loads((ROOT / 'download_manifest.json').read_text(encoding='utf-8')), key=lambda x: x['sample'])
    assert len(records) == 6 and pd.Series([r['condition'] for r in records]).value_counts().to_dict() == {'MCAO': 3, 'Sham': 3}
    audits = [aggregate(r) for r in records]
    m = pd.DataFrame([{'sample': r['sample'], 'condition': r['condition']} for r in records])
    allcounts = pd.concat([pd.read_csv(FOLDER / (r['sample'] + '_sample_counts.tsv'), sep='\t') for r in records], ignore_index=True)
    allcounts.to_csv(FOLDER / 'sample_cell_counts.tsv', sep='\t', index=False)
    allcounts.to_csv(OUT / 'third_cell_type_library_counts.tsv', sep='\t', index=False)
    mean, frac = {}, {}
    gene_sets = []
    coverage = []
    for r in records:
        s = r['sample']
        mean[s] = pd.read_csv(FOLDER / (s + '_mean_logcp10k.tsv.gz'), sep='\t', index_col=0)
        frac[s] = pd.read_csv(FOLDER / (s + '_fractions.tsv.gz'), sep='\t', index_col=0)
        assert mean[s].index.equals(frac[s].index)
        gene_sets.append(set(mean[s].index))
        coverage.append(pd.Series(True, index=mean[s].index, name=s))
    common_genes = set.intersection(*gene_sets)
    pd.concat(coverage, axis=1).fillna(False).to_csv(OUT / 'third_gene_feature_coverage.tsv.gz', sep='\t', index_label='gene')
    lrpath = BASE / 'literature' / 'mouseconsensus.csv'
    lr = pd.read_csv(lrpath).drop_duplicates(['source_genesymbol', 'target_genesymbol'])
    case = (m.condition == 'MCAO').to_numpy()
    samples = m['sample'].tolist()
    edge_rows, sample_rows = [], []
    for sender, receiver in [('Astrocyte', 'Endothelial'), ('Endothelial', 'Astrocyte')]:
        for ligand, receptor in lr[['source_genesymbol', 'target_genesymbol']].itertuples(index=False, name=None):
            lg, rg = ligand.split('_'), receptor.split('_')
            if not set(lg + rg) <= common_genes:
                continue
            l = np.array([mean[s].loc[lg, 'primary__' + sender].min() for s in samples])
            r = np.array([mean[s].loc[rg, 'primary__' + receiver].min() for s in samples])
            ld = np.array([frac[s].loc[lg, 'primary__' + sender].min() for s in samples])
            rd = np.array([frac[s].loc[rg, 'primary__' + receiver].min() for s in samples])
            score = np.sqrt(l * r)
            rec = dict(zip(KEY, (sender, receiver, ligand, receptor)))
            rec.update(dataset=ACC, config='primary', mean_control_score=float(score[~case].mean()),
                       mean_stroke_score=float(score[case].mean()), score_difference=float(score[case].mean() - score[~case].mean()))
            for threshold in (.05, .1, .2):
                tag = f'{int(threshold * 100):02d}'
                rec['eligible_' + tag] = bool(max((ld[case] >= threshold).sum(), (ld[~case] >= threshold).sum()) >= 2 and max((rd[case] >= threshold).sum(), (rd[~case] >= threshold).sum()) >= 2)
            edge_rows.append(rec)
            for j, s in enumerate(samples):
                sample_rows.append({**dict(zip(KEY, (sender, receiver, ligand, receptor))), 'dataset': ACC, 'sample': s,
                                   'condition': m.loc[j, 'condition'], 'score': score[j], 'ligand_logexpr': l[j],
                                   'receptor_logexpr': r[j], 'ligand_detection': ld[j], 'receptor_detection': rd[j],
                                   'sample_both_detected_10': bool(ld[j] >= .1 and rd[j] >= .1)})
    edges = pd.DataFrame(edge_rows)
    per_sample = pd.DataFrame(sample_rows)
    assert not edges.duplicated(KEY).any()
    edges.to_csv(OUT / 'third_all_covered_coavailability_effects.tsv.gz', sep='\t', index=False)
    per_sample.to_csv(OUT / 'third_all_covered_coavailability_library_scores.tsv.gz', sep='\t', index=False)
    summaries = compare(edges)
    ex = pd.DataFrame(EXAMPLES, columns=KEY).merge(edges, on=KEY, how='left', validate='1:1')
    ex.to_csv(OUT / 'third_fixed_examples_effects_eligibility.tsv', sep='\t', index=False)
    pd.DataFrame(EXAMPLES, columns=KEY).merge(per_sample, on=KEY, how='left', validate='1:m').to_csv(OUT / 'third_fixed_examples_library_scores.tsv', sep='\t', index=False)
    components = []
    for sender, receiver, ligand, receptor in EXAMPLES:
        for role, cell_type, genes in [('ligand', sender, ligand.split('_')), ('receptor', receiver, receptor.split('_'))]:
            for gene in genes:
                for j, s in enumerate(samples):
                    present = gene in mean[s].index
                    components.append({'sender': sender, 'receiver': receiver, 'ligand': ligand, 'receptor': receptor,
                                       'component_role': role, 'cell_type': cell_type, 'gene': gene, 'sample': s,
                                       'condition': m.loc[j, 'condition'], 'gene_in_submitted_features': present,
                                       'mean_logcp10k': float(mean[s].loc[gene, 'primary__' + cell_type]) if present else np.nan,
                                       'detection_fraction': float(frac[s].loc[gene, 'primary__' + cell_type]) if present else np.nan})
    pd.DataFrame(components).to_csv(OUT / 'third_fixed_examples_component_expression.tsv', sep='\t', index=False)
    metadata = {'dataset': ACC, 'n_libraries': 6, 'n_verified_independent_animals': None,
                'author_reported_biological_replicates_per_condition': 3,
                'common_submitted_symbols_all_six': len(common_genes), 'all_gene_covered_edges': len(edges),
                'eligible_05': int(edges.eligible_05.sum()), 'eligible_10': int(edges.eligible_10.sum()), 'eligible_20': int(edges.eligible_20.sum()),
                'formula': 'sqrt(min constituent mean log1p(CP10k) ligand * min constituent mean log1p(CP10k) receptor); delta is equal-library MCAO mean minus Sham mean.',
                'no_silent_missing_gene_imputation': True, 'cohort_p_values_computed': False,
                'does_not_enter_original_200_label_null': True,
                'resource_sha256': hashlib.sha256(lrpath.read_bytes()).hexdigest(), 'aggregation_checks': audits,
                'comparisons': summaries}
    (OUT / 'third_coavailability_audit.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(metadata, ensure_ascii=False, indent=2), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary-only', action='store_true')
    args = parser.parse_args()
    if args.summary_only:
        summaries = compare(pd.read_csv(OUT/'third_all_covered_coavailability_effects.tsv.gz', sep='\t'))
        path = OUT/'third_coavailability_audit.json'
        audit = json.loads(path.read_text(encoding='utf-8'))
        audit['comparisons'] = summaries
        audit['author_reported_biological_replicates_per_condition'] = 3
        path.write_text(json.dumps(audit,indent=2),encoding='utf-8')
        print(json.dumps(summaries,indent=2))
    else:
        main()
