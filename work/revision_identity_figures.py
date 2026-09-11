"""Two identity-review figures; run only after all 11 final sample audits exist.

No cell labels, thresholds, or original figures are changed. Figure 6 uses
animal-level proportions. The supplementary histograms are technical QC only.
"""
from pathlib import Path
import argparse
import json
import re
import string

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'revision_results' / 'identity'
OUT = ROOT.parent / 'outputs' / 'revised' / 'figures'
DATA = ROOT.parent / 'outputs' / 'revised' / 'figure_source_data'
COHORTS = ('GSE174574', 'GSE245386')
EXPECTED = {
    'GSE174574': [f'GSM{i}' for i in range(5319987, 5319993)],
    'GSE245386': [f'GSM{i}' for i in range(7841720, 7841725)],
}
COHORT_NAMES = {'GSE174574': 'Discovery', 'GSE245386': 'Validation'}
BLUE, ORANGE, GRAY = '#0072B2', '#D55E00', '#8A8A8A'
TARGETS = ('Astrocyte', 'Endothelial')
CONFIDENCE = 0.5  # Same rule as revision_aggregate.py; not a new decision rule.
BUCKETS = ('Astrocyte', 'Endothelial', 'Other assigned', 'Unassigned')
METRICS = (
    ('whole_brain_agreement', 'Whole-brain match', 'o', -0.20, True),
    ('cortex_hippocampus_agreement', 'Cortex/hippocampus match', 'D', 0.00, True),
    ('reference_singlet_retained', 'High-confidence matched singlets', 's', 0.20, False),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def boolean_column(series, name):
    values = series.astype(str).str.lower()
    require(values.isin(['true', 'false', '1', '0']).all(), f'Invalid boolean: {name}')
    return values.isin(['true', '1']).to_numpy()


def sample_metadata():
    """Check frozen sample IDs against original QC and submitted file names."""
    records = []
    manifest = pd.read_csv(ROOT / 'revision_sample_manifest.tsv', sep='\t').set_index(['dataset', 'sample'])
    require(manifest.index.is_unique, 'Duplicate frozen manifest sample')
    for acc in COHORTS:
        qc = pd.read_csv(ROOT / 'processed' / acc / 'qc_summary.tsv', sep='\t')
        require(qc['sample'].is_unique, f'Duplicate QC sample IDs: {acc}')
        require(set(qc['sample']) == set(EXPECTED[acc]), f'Unexpected QC samples: {acc}')
        qc = qc.set_index('sample')
        for sample in EXPECTED[acc]:
            matrices = list((ROOT / 'data' / acc / 'raw').glob(f'{sample}_*matrix.mtx.gz'))
            prefix = manifest.loc[(acc, sample), 'source_prefix']
            if matrices:
                require(len(matrices) == 1 and matrices[0].name.removesuffix('_matrix.mtx.gz') == prefix,
                        f'Original filename disagrees with frozen manifest: {sample}')
            suffix = prefix[len(sample) + 1:]
            match = re.fullmatch(r'(sham|MCAO|WTC|WTM)(\d+)', suffix, flags=re.I)
            require(match is not None, f'Unrecognized condition suffix: {prefix}')
            condition = 'Sham' if match.group(1).lower() in ('sham', 'wtc') else 'MCAO'
            records.append({
                'dataset': acc, 'sample': sample, 'condition': condition,
                'sample_alias': f'{condition}{match.group(2)}', 'source_prefix': prefix,
                'n_qc': int(qc.loc[sample, 'n_qc']),
                'n_rule_astrocyte': int(qc.loc[sample, 'Astrocyte']),
                'n_rule_endothelial': int(qc.loc[sample, 'Endothelial']),
            })
    meta = pd.DataFrame(records)
    observed = meta.groupby(['dataset', 'condition']).size().to_dict()
    require(observed == {('GSE174574', 'Sham'): 3, ('GSE174574', 'MCAO'): 3,
                         ('GSE245386', 'Sham'): 3, ('GSE245386', 'MCAO'): 2},
            'Condition counts do not match the study design')
    return meta


def load_complete_inputs():
    """Validate every input before creating any user-facing output."""
    meta = sample_metadata()
    missing = []
    for row in meta.itertuples(index=False):
        for suffix in ('_audit.json', '_identity.tsv.gz', '_simulated_doublet_scores.tsv.gz'):
            path = INPUT / row.dataset / (row.sample + suffix)
            if not path.is_file():
                missing.append(str(path))
    if missing:
        raise FileNotFoundError('All 11 final audits and score files are required. Missing:\n'
                                + '\n'.join(missing))
    frames, audits, simulated = {}, {}, {}
    required = ['barcode', 'sample', 'cell_type', 'whole_brain_broad',
                'whole_brain_probability', 'cortex_hippocampus_broad',
                'cortex_hippocampus_probability', 'predicted_doublet', 'doublet_score']
    for row in meta.itertuples(index=False):
        key = (row.dataset, row.sample)
        stem = INPUT / row.dataset / row.sample
        audit = json.loads(Path(str(stem) + '_audit.json').read_text(encoding='utf-8-sig'))
        frame = pd.read_csv(str(stem) + '_identity.tsv.gz', sep='\t')
        sim = pd.read_csv(str(stem) + '_simulated_doublet_scores.tsv.gz', sep='\t')
        require(set(required).issubset(frame.columns), f'Missing identity columns: {key}')
        require(frame['barcode'].is_unique and frame['barcode'].notna().all(),
                f'Invalid barcodes: {key}')
        require(frame['sample'].eq(row.sample).all(), f'Sample mismatch: {key}')
        require(audit['dataset'] == row.dataset and audit['sample'] == row.sample,
                f'Audit identity mismatch: {key}')
        require(len(frame) == row.n_qc == audit['n_qc'], f'QC count mismatch: {key}')
        for field in ('cell_type', 'whole_brain_broad', 'cortex_hippocampus_broad'):
            require(frame[field].notna().all(), f'Missing label in {field}: {key}')
        for field in ('whole_brain_probability', 'cortex_hippocampus_probability', 'doublet_score'):
            values = frame[field].to_numpy(dtype=float)
            require(np.isfinite(values).all() and ((values >= 0) & (values <= 1)).all(),
                    f'Invalid score in {field}: {key}')
        frame['predicted_doublet'] = boolean_column(frame['predicted_doublet'], str(key))
        scrub = audit['scrublet']
        threshold = float(scrub['threshold'])
        require(np.isfinite(threshold) and 0 <= threshold <= 1, f'Invalid threshold: {key}')
        require(scrub.get('automatic_threshold') is True, f'Expected automatic threshold: {key}')
        require(int(frame['predicted_doublet'].sum()) == scrub['n_predicted'],
                f'Doublet count mismatch: {key}')
        require(np.array_equal(frame['doublet_score'].to_numpy() > threshold,
                               frame['predicted_doublet'].to_numpy()),
                f'Predictions differ from stored threshold: {key}')
        require('simulated_score' in sim, f'Missing simulated scores: {key}')
        values = sim['simulated_score'].to_numpy(dtype=float)
        require(len(values) == scrub['n_simulated'] and len(values) > 0,
                f'Simulated count mismatch: {key}')
        require(np.isfinite(values).all() and ((values >= 0) & (values <= 1)).all(),
                f'Invalid simulated scores: {key}')
        for typ, count in zip(TARGETS, (row.n_rule_astrocyte, row.n_rule_endothelial)):
            require(int(frame['cell_type'].eq(typ).sum()) == count, f'Target QC mismatch: {key}, {typ}')
            require(count > 0, f'Zero denominator for original target cells: {key}, {typ}')
        frames[key], audits[key], simulated[key] = frame, audit, values
    return meta, frames, audits, simulated


def configure_style():
    plt.rcParams.update({
        'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9,
        'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'axes.spines.top': False,
        'axes.spines.right': False, 'savefig.dpi': 400,
    })


def save(fig, name):
    for extension in ('png', 'pdf', 'svg'):
        fig.savefig(OUT / f'{name}.{extension}', bbox_inches='tight', dpi=400)
    plt.close(fig)
    print(name, flush=True)


def panel_letter(ax, letter):
    ax.text(-0.13, 1.06, letter, transform=ax.transAxes, weight='bold', fontsize=10)


def coarse(series):
    return series.where(series.isin(('Astrocyte', 'Endothelial', 'Unassigned')), 'Other assigned')


def identity_tables(meta, frames):
    proportions, full_confusion, coarse_confusion = [], [], []
    for row in meta.itertuples(index=False):
        frame = frames[(row.dataset, row.sample)]
        base = {'dataset': row.dataset, 'sample': row.sample, 'condition': row.condition,
                'sample_alias': row.sample_alias}
        for typ in TARGETS:
            original = frame['cell_type'].eq(typ)
            n = int(original.sum())
            whole_match = frame['whole_brain_broad'].eq(typ)
            retained = (whole_match & frame['whole_brain_probability'].ge(CONFIDENCE)
                        & ~frame['predicted_doublet'])
            for metric, mask in (
                ('whole_brain_agreement', whole_match),
                ('cortex_hippocampus_agreement', frame['cortex_hippocampus_broad'].eq(typ)),
                ('reference_singlet_retained', retained),
            ):
                numerator = int((original & mask).sum())
                proportions.append({**base, 'original_cell_type': typ, 'metric': metric,
                                    'numerator': numerator, 'denominator': n,
                                    'fraction': numerator / n})
        for reference in ('whole_brain', 'cortex_hippocampus'):
            full = frame.groupby(['cell_type', reference + '_broad'], observed=True).size()
            for (original, predicted), count in full.items():
                full_confusion.append({**base, 'reference': reference,
                                       'original_label': original, 'reference_label': predicted,
                                       'n_cells': int(count)})
        table = pd.crosstab(coarse(frame['cell_type']), coarse(frame['whole_brain_broad']))
        table = table.reindex(index=BUCKETS, columns=BUCKETS, fill_value=0)
        require(int(table.to_numpy().sum()) == len(frame), 'Coarse confusion dropped cells')
        for original in BUCKETS:
            denom = int(table.loc[original].sum())
            for predicted in BUCKETS:
                count = int(table.loc[original, predicted])
                coarse_confusion.append({**base, 'original_bucket': original,
                                         'reference_bucket': predicted, 'n_cells': count,
                                         'row_denominator': denom,
                                         'row_fraction': count / denom if denom else np.nan})
    return pd.DataFrame(proportions), pd.DataFrame(full_confusion), pd.DataFrame(coarse_confusion)


def figure6(meta, frame, confusion):
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 8.1), layout='constrained',
                             gridspec_kw={'height_ratios': [1, 1, 1.1]})
    for row_index, typ in enumerate(TARGETS):
        for col_index, acc in enumerate(COHORTS):
            ax = axes[row_index, col_index]
            rows = meta[meta['dataset'].eq(acc)]
            subset = frame[frame['dataset'].eq(acc) & frame['original_cell_type'].eq(typ)]
            for x, animal in enumerate(rows.itertuples(index=False)):
                color = BLUE if animal.condition == 'Sham' else ORANGE
                for metric, _, marker, offset, filled in METRICS:
                    value = subset.loc[subset['sample'].eq(animal.sample) & subset['metric'].eq(metric), 'fraction']
                    require(len(value) == 1, 'Each animal/metric needs exactly one value')
                    ax.scatter(x + offset, 100 * value.iloc[0], marker=marker, s=27,
                               facecolors=color if filled else 'none', edgecolors=color,
                               linewidths=0.85, zorder=3)
            ax.set_xticks(range(len(rows)), rows['sample_alias'], rotation=35, ha='right')
            ax.set_ylim(0, 103)
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.set_xlim(-0.55, len(rows) - 0.45)
            ax.set_ylabel('Original target cells (%)')
            ax.set_title(f'{COHORT_NAMES[acc]}: {typ}')
            ax.grid(axis='y', color='#E5E5E5', linewidth=0.5)
            ax.set_axisbelow(True)
            panel_letter(ax, string.ascii_uppercase[2 * row_index + col_index])
    means = confusion.groupby(['dataset', 'original_bucket', 'reference_bucket'], observed=True)['row_fraction'].agg(['mean', 'count']).reset_index()
    ticklabels = ['Astrocyte', 'Endothelial', 'Other assigned', 'Unassigned']
    for index, acc in enumerate(COHORTS):
        ax = axes[2, index]
        block = means[means['dataset'].eq(acc)]
        matrix = block.pivot(index='original_bucket', columns='reference_bucket', values='mean').reindex(index=BUCKETS, columns=BUCKETS)
        im = ax.imshow(matrix, vmin=0, vmax=1, cmap='cividis', aspect='auto')
        ax.set_xticks(range(4), ticklabels, rotation=40, ha='right')
        ax.set_yticks(range(4), ticklabels)
        ax.set_xlabel('Whole-brain reference label')
        ax.set_ylabel('Original rule label')
        ax.set_title(f'{COHORT_NAMES[acc]}: all QC cells')
        for i in range(4):
            for j in range(4):
                value = matrix.iloc[i, j]
                ax.text(j, i, f'{100 * value:.0f}' if np.isfinite(value) else 'NA',
                        ha='center', va='center', fontsize=7,
                        color='black' if np.isfinite(value) and value > .55 else 'white')
        panel_letter(ax, string.ascii_uppercase[4 + index])
    fig.colorbar(im, ax=axes[2, :], shrink=.65, label='Mean within-animal row fraction')
    metric_handles = [Line2D([0], [0], marker=m, linestyle='none', color='black',
                             markerfacecolor='black' if fill else 'none', markersize=4,
                             label=title) for _, title, m, _, fill in METRICS]
    fig.legend(handles=metric_handles, loc='outside lower center', ncol=1, frameon=False)
    means.to_csv(DATA / 'figure6_confusion_animal_mean.tsv', sep='\t', index=False)
    save(fig, 'Figure6_independent_identity_review')


def figure_s3(meta, frames, audits, simulated):
    # Common bin edges include every finite score and every stored threshold.
    maximum = max(max(float(frames[key]['doublet_score'].max()), float(values.max()),
                      float(audits[key]['scrublet']['threshold']))
                  for key, values in simulated.items())
    upper = min(1., max(.1, np.ceil(maximum * 20) / 20))
    bins = np.linspace(0., upper, 61)
    fig, axes = plt.subplots(4, 3, figsize=(7.2, 8.4), layout='constrained', sharex=True)
    bin_records, audit_records = [], []
    for k, row in enumerate(meta.itertuples(index=False)):
        key = (row.dataset, row.sample)
        ax = axes.flat[k]
        observed = frames[key]['doublet_score'].to_numpy()
        scrub = audits[key]['scrublet']
        for name, values, color, linestyle in (
            ('Observed', observed, BLUE, '-'),
            ('Simulated', simulated[key], ORANGE, '--'),
        ):
            counts, _ = np.histogram(values, bins=bins)
            require(int(counts.sum()) == len(values), f'Histogram range dropped scores: {key}')
            fractions = counts / len(values)
            ax.stairs(fractions, bins, color=color, linestyle=linestyle, linewidth=.9)
            for left, right, count, fraction in zip(bins[:-1], bins[1:], counts, fractions):
                bin_records.append({'dataset': row.dataset, 'sample': row.sample,
                                    'condition': row.condition, 'distribution': name,
                                    'bin_left': left, 'bin_right': right, 'n_scores': int(count),
                                    'total_scores': len(values), 'fraction_per_bin': fraction})
        threshold = float(scrub['threshold'])
        ax.axvline(threshold, color='black', linestyle='-.', linewidth=.8)
        ax.set_xlim(0, upper)
        ax.set_ylim(bottom=0)
        ax.set_title(f'{COHORT_NAMES[row.dataset]}: {row.sample_alias}\n{row.sample}; n={len(observed):,}', fontsize=8)
        ax.text(.97, .96, f't={threshold:.3f}\nflagged {100 * scrub["predicted_fraction"]:.1f}%',
                transform=ax.transAxes, ha='right', va='top', fontsize=6.5,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .8, 'pad': 1})
        ax.tick_params(axis='x', labelbottom=True)
        ax.set_xlabel('Scrublet score')
        if k % 3 == 0:
            ax.set_ylabel('Fraction per bin')
        panel_letter(ax, string.ascii_uppercase[k])
        audit_records.append({'dataset': row.dataset, 'sample': row.sample,
                              'condition': row.condition, 'n_observed': len(observed), **scrub})
    axes.flat[-1].set_axis_off()
    handles = [Line2D([0], [0], color=BLUE, lw=1, label='Observed cells'),
               Line2D([0], [0], color=ORANGE, lw=1, ls='--', label='Simulated doublets'),
               Line2D([0], [0], color='black', lw=1, ls='-.', label='Stored automatic threshold')]
    axes.flat[-1].legend(handles=handles, loc='center', frameon=False, fontsize=7)
    axes.flat[-1].text(.5, .18, 'Per-library technical QC\nNo cell-level group P values',
                       transform=axes.flat[-1].transAxes, ha='center', va='center', fontsize=7)
    pd.DataFrame(bin_records).to_csv(DATA / 'figureS3_histogram_bins.tsv', sep='\t', index=False)
    pd.DataFrame(audit_records).to_csv(DATA / 'figureS3_sample_doublet_audits.tsv', sep='\t', index=False)
    save(fig, 'FigureS3_doublet_score_distributions')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-inputs', action='store_true', help='Validate all 11 inputs without generating figures')
    args = parser.parse_args()
    meta, frames, audits, simulated = load_complete_inputs()
    if args.check_inputs:
        print(f'All {len(meta)} sample audits and score matrices are complete and consistent.')
        return
    proportions, full, coarse_table = identity_tables(meta, frames)
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    meta.to_csv(DATA / 'figure6_sample_metadata.tsv', sep='\t', index=False)
    proportions.to_csv(DATA / 'figure6_target_proportions.tsv', sep='\t', index=False)
    full.to_csv(DATA / 'figure6_full_label_confusion.tsv', sep='\t', index=False)
    coarse_table.to_csv(DATA / 'figure6_coarse_confusion_by_sample.tsv', sep='\t', index=False)
    figure6(meta, proportions, coarse_table)
    figure_s3(meta, frames, audits, simulated)
    captions = """Figure 6. Independent reference-label review and retention of original target cells.
A-D: Each animal has three points. Filled circles/diamonds indicate the percentage of the original rule-defined target cells assigned the same broad label by the whole-brain/cortex-hippocampus reference, irrespective of confidence. Open squares indicate the fraction that also has whole-brain confidence >=0.5 and is not flagged by Scrublet, matching the reference_singlet sensitivity definition. All denominators are the original target cells in that animal; full numerators, denominators, sample aliases and GSM IDs are supplied. Blue indicates Sham and orange indicates MCAO. Percentages are descriptive and are not classification accuracy against ground truth.
E-F: All QC cells, including original Unassigned cells and every non-target label, enter the label comparison. Other assigned pools the remaining labels for display; the full original/reference label table is provided separately. Heatmap values are within-animal row percentages averaged equally across animals with a nonzero row denominator; the numbers in cells are percentages. Low-confidence best-match assignments are retained in these matrices. No group significance tests are performed. Reference confidence is a model score, and healthy-reference bias remains possible.

Figure S3. Per-library observed and simulated Scrublet score distributions.
Solid blue and dashed orange step histograms show observed cells and simulated doublets, respectively; colors in this figure denote distribution type. A black dash-dot line marks the stored automatic threshold without manual changes. All 11 libraries share bin edges and x-axis limits; each distribution is normalized by its own number of scores. Per-panel y axes may differ to keep distribution shapes visible. Titles identify cohort, treatment, GSM ID and number of observed QC cells; annotations give the stored threshold and flagged fraction. These are technical QC distributions, not animal-level treatment comparisons; no cell-level group P values are computed. Doublet exclusion does not constitute ambient RNA correction.
"""
    (DATA / 'figure6_and_S3_captions.txt').write_text(captions, encoding='utf-8')
    metadata = {'n_libraries': len(meta), 'n_qc_cells': int(meta['n_qc'].sum()),
                'confidence_threshold': CONFIDENCE, 'dpi': 400, 'font': 'Arial',
                'pdf_fonttype': 42, 'svg_fonttype': 'none',
                'figure6_analysis_unit': 'animal for proportions and equal-weight heatmap means',
                'doublet_thresholds_changed': False, 'cell_labels_changed': False,
                'ambient_RNA_corrected_by_this_script': False,
                'visual_review_status': 'Requires visual inspection after rendering'}
    (DATA / 'identity_figure_metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
