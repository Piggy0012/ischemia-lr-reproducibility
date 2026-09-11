"""Figure 7: actual cross-cohort selection/method robustness, no imputed data.

First run writes the PNG preview and machine QA. After inspecting the preview,
rerun with --export to create the final PNG/PDF/SVG at exactly 7.2 x 8.15 inches.
Incomplete summary inputs raise before any figure/source-data is written.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
from PIL import Image

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / 'outputs/revised'
NAME = 'Figure7_method_and_selection_robustness'
BLUE, ORANGE, GRAY = '#0072B2', '#D55E00', '#8A8A8A'
CONFIGS = ['primary', 'singlet', 'reference_singlet', 'reference_only']
CONFIG_LABELS = ['Original rules', 'Rules + singlets',
                 'Reference-supported singlets', 'Reference-only singlets']
LIANA_CONFIGS = ['primary', 'reference_singlet']
METRICS = ['lr_means', 'expr_prod', 'lr_logfc', 'spec_weight', 'lrscore', 'magnitude_rank']
METRIC_LABELS = ['CellPhoneDB mean', 'Expression product', 'Cell-type logFC',
                 'NATMI specificity', 'SingleCellSignalR', 'LIANA 1 − magnitude rank']
KEY = ['sender', 'receiver', 'ligand', 'receptor']
EXAMPLES = [
    ('Astrocyte', 'Endothelial', 'Timp3', 'Kdr'),
    ('Endothelial', 'Astrocyte', 'Ptn', 'Ptprz1'),
    ('Endothelial', 'Astrocyte', 'Plat', 'Lrp1'),
    ('Astrocyte', 'Endothelial', 'Spp1', 'Itga5_Itgb1'),
    ('Astrocyte', 'Endothelial', 'Col4a1', 'Itga3_Itgb1'),
]
DATASETS = ['GSE174574', 'GSE245386']
SKILL = Path('C:/Users/admin/.codex/skills/scipilot-figure-skill/scripts')


def bool_value(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if str(value).lower() in ['true', '1']:
        return True
    if str(value).lower() in ['false', '0'] or pd.isna(value):
        return False
    raise ValueError(f'Invalid eligibility boolean: {value!r}')


def safe_json(value):
    if isinstance(value, dict):
        return {str(k): safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [safe_json(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def load_inputs(table_dir):
    paths = {n: table_dir / n for n in ['summary.json', 'custom_cross_cohort_summary.tsv',
        'liana_cross_cohort_summary.tsv', 'custom_fixed_examples.tsv',
        'liana_fixed_examples_animal_effects.tsv']}
    absent = [str(p) for p in paths.values() if not p.exists()]
    if absent:
        raise FileNotFoundError('Figure 7 awaits completed analysis tables; no figure generated:\n' + '\n'.join(absent))
    summary = json.loads(paths['summary.json'].read_text(encoding='utf-8'))
    if summary.get('status') != 'complete':
        raise ValueError('Figure 7 requires summary.json status=complete; no partial or fabricated figure is allowed')
    a = pd.read_csv(paths['custom_cross_cohort_summary.tsv'], sep='\t')
    a = a[np.isclose(a.expression_gate, .1)].set_index('config').loc[CONFIGS].reset_index()
    b = pd.read_csv(paths['liana_cross_cohort_summary.tsv'], sep='\t')
    assert not b.duplicated(['config', 'metric']).any()
    b = b.set_index(['config', 'metric']).loc[pd.MultiIndex.from_product([LIANA_CONFIGS, METRICS])].reset_index()
    for d in [a, b]:
        assert (d.n_same_direction >= 0).all()
        assert (d.n_same_direction <= d.n_direction_evaluable).all()
        assert (d.n_direction_evaluable <= d.n_common).all()
        valid = d.n_direction_evaluable > 0
        assert np.allclose(d.loc[valid, 'same_direction_fraction'],
                           d.loc[valid, 'n_same_direction'] / d.loc[valid, 'n_direction_evaluable'])
        assert d.loc[valid, 'same_direction_fraction'].between(0, 1).all()
    custom = pd.read_csv(paths['custom_fixed_examples.tsv'], sep='\t').set_index(['dataset', 'config'] + KEY)
    liana = pd.read_csv(paths['liana_fixed_examples_animal_effects.tsv'], sep='\t')
    liana = liana[liana.metric == 'magnitude_rank'].set_index(['dataset', 'config'] + KEY)
    assert custom.index.is_unique and liana.index.is_unique
    matrix_rows = []
    columns = [('custom', c) for c in CONFIGS] + [('liana', c) for c in LIANA_CONFIGS]
    for r, example in enumerate(EXAMPLES):
        for col, (method, config) in enumerate(columns):
            for acc in DATASETS:
                index = (acc, config) + example
                d = custom if method == 'custom' else liana
                if index not in d.index:
                    raise ValueError(f'Fixed example absent from completed summary: {index}, {method}')
                row = d.loc[index]
                value = row.get('score_difference', np.nan)
                if method == 'custom':
                    eligible = bool_value(row.get('eligible_10', False))
                else:
                    eligible = str(row.status).startswith('complete_case')
                evaluable = eligible and np.isfinite(value)
                sign = '+' if evaluable and value > 0 else '−' if evaluable and value < 0 else '0' if evaluable else 'NA'
                matrix_rows.append(dict(zip(KEY, example)) | {
                    'row': r, 'column': col, 'method': method, 'config': config, 'dataset': acc,
                    'source_status': row.status, 'evaluable': evaluable,
                    'source_score_difference': value, 'display_sign': sign,
                    'custom_eligible_10': row.get('eligible_10', np.nan),
                    'liana_missing_animals': row.get('missing_or_nonfinite_samples', '')})
    c = pd.DataFrame(matrix_rows)
    assert len(a) == 4 and len(b) == 12 and len(c) == 60
    source_hashes = {n: hashlib.sha256(p.read_bytes()).hexdigest() for n, p in paths.items()}
    return a, b, c, source_hashes


def profile_inputs(a, b, c, qa_dir):
    # Candidate summaries are descriptive proportions, not animal observations.
    profile = {'purpose': 'Describe sensitivity of cross-cohort signs to cell selection and scoring method',
        'unit_AB': 'one statistic per evaluated candidate set; candidates share ligand/receptor components',
        'unit_C': 'one fixed candidate / method-selection / cohort effect sign',
        'A_rows': len(a), 'B_rows': len(b), 'C_rows': len(c),
        'A_candidate_denominators': a[['config', 'n_common', 'n_direction_evaluable']].to_dict('records'),
        'B_candidate_denominators': b[['config', 'metric', 'n_common', 'n_direction_evaluable']].to_dict('records'),
        'C_missing_by_method': c.groupby('method').evaluable.agg(['size', 'sum']).to_dict(),
        'chart_choice': 'Aligned dot plots with explicit n/d; categorical sign-and-eligibility matrix',
        'alternatives_considered': ['Annotated numeric table', 'Separate sign matrices by cohort'],
        'not_used': 'No binomial CI, correlation significance, pooled-cell error bars, or mean bars',
        'cross_cohort_n_animals': {'discovery': {'Sham': 3, 'MCAO': 3}, 'validation': {'Sham': 3, 'MCAO': 2}}}
    if SKILL.exists():
        sys.path.insert(0, str(SKILL))
        from profile_data import profile_data
        profile['scipilot_A'] = profile_data(a, group_cols=['config'])
        profile['scipilot_B'] = profile_data(b, group_cols=['config', 'metric'])
    qa_dir.mkdir(parents=True, exist_ok=True)
    (qa_dir / 'data_profile.json').write_text(json.dumps(safe_json(profile), indent=2,
                                                       default=str, ensure_ascii=False), encoding='utf-8')


def style_dots(ax, n, labels):
    ax.set_xlim(0, 100)
    ax.set_ylim(-.6, n - .4)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks(np.arange(n), labels[::-1])
    ax.grid(axis='x', color='#E5E5E5', lw=.55, zorder=0)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=7)
    ax.set_xlabel('Cross-cohort same-direction candidates (%)', labelpad=4)


def make_figure(a, b, c):
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9,
        'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False,
        'savefig.dpi': 400})
    fig = plt.figure(figsize=(7.2, 8.15))
    fig.text(.035, .965, 'A', fontweight='bold', fontsize=10)
    fig.text(.083, .965, 'Cell-selection sensitivity of custom scores', fontsize=9)
    ax = fig.add_axes([.275, .785, .40, .14])
    ann = fig.add_axes([.715, .785, .26, .14], sharey=ax)
    ann.set_xlim(0, 1); ann.set_axis_off()
    style_dots(ax, 4, CONFIG_LABELS)
    ann.text(.02, 1.12, 'Same / evaluable', fontsize=7, transform=ann.transAxes)
    ann.text(.86, 1.12, 'rho', fontsize=7, transform=ann.transAxes, ha='center')
    for i, config in enumerate(CONFIGS):
        row = a[a.config == config].iloc[0]
        y = 3 - i
        if np.isfinite(row.same_direction_fraction):
            ax.scatter(row.same_direction_fraction * 100, y, color=BLUE, marker='o', s=28, zorder=3)
        else:
            ax.text(50, y, 'Not evaluable', ha='center', va='center', fontsize=7)
        ann.text(.02, y, f'{int(row.n_same_direction)}/{int(row.n_direction_evaluable)}', va='center', fontsize=7.5)
        rho = f'{row.spearman_rho:.2f}' if np.isfinite(row.spearman_rho) else 'NA'
        ann.text(.86, y, rho, va='center', ha='center', fontsize=7.5)

    fig.text(.035, .733, 'B', fontweight='bold', fontsize=10)
    fig.text(.083, .733, 'Method sensitivity within the same data', fontsize=9)
    handles = [Line2D([0], [0], marker='o', color=BLUE, linestyle='none', markersize=4, label='Original rules'),
               Line2D([0], [0], marker='^', color=ORANGE, linestyle='none', markersize=4, label='Reference-supported singlets')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(.265, .69), frameon=False,
               ncol=2, handletextpad=.35, columnspacing=1)
    ax = fig.add_axes([.275, .46, .40, .222])
    ann = fig.add_axes([.715, .46, .26, .222], sharey=ax)
    ann.set_xlim(0, 1); ann.set_axis_off()
    style_dots(ax, 6, METRIC_LABELS)
    ann.text(.02, 1.08, 'Same / evaluable', fontsize=7, transform=ann.transAxes)
    for i, metric in enumerate(METRICS):
        for j, config in enumerate(LIANA_CONFIGS):
            row = b[(b.config == config) & (b.metric == metric)].iloc[0]
            y = 5 - i + (.145 if j == 0 else -.145)
            color, marker = (BLUE, 'o') if j == 0 else (ORANGE, '^')
            if np.isfinite(row.same_direction_fraction):
                ax.scatter(row.same_direction_fraction * 100, y, color=color, marker=marker, s=23, zorder=3)
            else:
                ax.text(50, y, 'Not evaluable', ha='center', va='center', fontsize=6.5, color=color)
            ann.scatter(.03, y, color=color, marker=marker, s=12, clip_on=False)
            ann.text(.14, y, f'{int(row.n_same_direction)}/{int(row.n_direction_evaluable)}',
                     color=color, va='center', fontsize=7.1)
    fig.text(.275, .401, 'Fractions describe candidate sets; no edge-level confidence intervals.', fontsize=6.6)

    fig.text(.035, .363, 'C', fontweight='bold', fontsize=10)
    fig.text(.083, .363, 'Fixed-candidate direction and eligibility audit', fontsize=9)
    ax = fig.add_axes([.34, .135, .63, .158])
    ax.set_xlim(-.5, 5.5); ax.set_ylim(4.5, -.5)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    ax.set_xticks(range(6), ['Rules', 'Singlets', 'Ref.-supported\nsinglets', 'Ref.-only\nsinglets',
                            'Rules', 'Ref.-supported\nsinglets'], fontsize=6.8)
    ax.xaxis.tick_top(); ax.tick_params(axis='both', length=0, pad=5)
    labels = [f'{ligand}–{receptor.replace("_", "/")}\n' + ('A → E' if sender == 'Astrocyte' else 'E → A')
              for sender, receiver, ligand, receptor in EXAMPLES]
    ax.set_yticks(range(5), labels, fontsize=7)
    for i in range(5):
        for j in range(6):
            q = c[(c.row == i) & (c.column == j)].set_index('dataset')
            signs = q.loc[DATASETS, 'display_sign'].tolist()
            if 'NA' in signs:
                fill = '#E6E6E6'
            elif signs == ['+', '+']:
                fill = '#F7D7C7'
            elif signs == ['−', '−']:
                fill = '#CEE5F2'
            else:
                fill = '#FFFFFF'
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor=fill,
                                   edgecolor='#C1C1C1', linewidth=.55))
            ax.text(j, i, ' / '.join(signs), ha='center', va='center', fontsize=7.7)
    ax.axvline(3.5, color='black', linewidth=1.15)
    fig.text(.55, .337, 'Custom coavailability', ha='center', fontsize=7.8)
    fig.text(.865, .337, 'LIANA 1 − rank', ha='center', fontsize=7.8)
    fig.text(.083, .103, 'Entries: discovery / validation effect signs; NA = not evaluable.', fontsize=7)
    fig.text(.083, .081, 'Custom: 10% cohort gate. LIANA: complete observed scores in every animal.', fontsize=7)
    fig.text(.083, .059, 'A → E: astrocyte to endothelium; E → A: endothelium to astrocyte.', fontsize=7)
    fig.legend(handles=[Patch(facecolor='#CEE5F2', edgecolor='#C1C1C1', label='Both decrease'),
        Patch(facecolor='#F7D7C7', edgecolor='#C1C1C1', label='Both increase'),
        Patch(facecolor='white', edgecolor='#C1C1C1', label='Discordant / zero'),
        Patch(facecolor='#E6E6E6', edgecolor='#C1C1C1', label='At least one NA')],
        loc='lower center', bbox_to_anchor=(.53, .012), frameon=False, ncol=4,
        handlelength=1.1, columnspacing=1.1, handletextpad=.45)
    return fig


def main(table_dir, output_base, export=False):
    a, b, c, hashes = load_inputs(table_dir)
    qa = ROOT / 'revision_figure7_qa'
    profile_inputs(a, b, c, qa)
    fig = make_figure(a, b, c)
    qa.mkdir(parents=True, exist_ok=True)
    issues = []
    if SKILL.exists():
        sys.path.insert(0, str(SKILL))
        from visual_qa import audit_layout
        issues = audit_layout(fig)
    fig.savefig(qa / 'preview.png', dpi=180, facecolor='white')
    with Image.open(qa / 'preview.png') as im:
        im.convert('L').save(qa / 'preview_grayscale.png')
    report = {'figure': NAME, 'inches': [7.2, 8.15], 'dpi': 400, 'layout_issues': issues,
              'input_sha256': hashes, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'n_panels': 3, 'human_visual_review_required': True, 'export_requested': export}
    (qa / 'qa.json').write_text(json.dumps(safe_json(report), indent=2, ensure_ascii=False), encoding='utf-8')
    if any(severity == 'FAIL' for severity, _ in issues):
        raise RuntimeError(f'Figure layout failed; inspect {qa / "qa.json"}')
    data = output_base / 'figure_source_data'
    data.mkdir(parents=True, exist_ok=True)
    a.to_csv(data / 'figure7_A_custom_concordance.tsv', sep='\t', index=False)
    b.to_csv(data / 'figure7_B_liana_concordance.tsv', sep='\t', index=False)
    c.to_csv(data / 'figure7_C_direction_eligibility.tsv', sep='\t', index=False)
    (data / 'figure7_caption.md').write_text(CAPTION, encoding='utf-8')
    if export:
        out = output_base / 'figures'
        out.mkdir(parents=True, exist_ok=True)
        for ext in ['png', 'pdf', 'svg']:
            fig.savefig(out / f'{NAME}.{ext}', dpi=400, facecolor='white')
        with Image.open(out / f'{NAME}.png') as im:
            assert im.size == (2880, 3260)
            assert min(im.info['dpi']) >= 399.9
        if SKILL.exists():
            from check_figure import check_figure
            report['png_file_check'] = check_figure(str(out / f'{NAME}.png'), min_dpi=400,
                                                    target_inches=(7.2, 8.15))
            (qa / 'qa.json').write_text(json.dumps(safe_json(report), indent=2, ensure_ascii=False,
                                                 default=str), encoding='utf-8')
    plt.close(fig)
    print(json.dumps({'preview': str(qa / 'preview.png'), 'layout_issues': issues,
                      'exported': export, 'source_data': str(data)}, ensure_ascii=False), flush=True)


CAPTION = '''**Figure 7. Sensitivity of cross-cohort candidate concordance to cell selection and communication methods.**

A, The four cell selections are original rule labels, rule labels excluding predicted doublets, reference-supported rule labels excluding predicted doublets, and reference-only target labels excluding predicted doublets. Points show same-sign fractions among candidates meeting the cohort-level 10% expression gate in both studies. Numerators, direction-evaluable denominators and descriptive Spearman correlations are printed alongside the points.

B, The same-sign fractions are calculated from actual per-animal LIANA scores, retaining only candidates with a finite expression-eligible score in every animal within each study and present in both study-level complete-case sets. Blue circles denote original rules; orange triangles denote reference-supported singlets. The expression-product column is shared by Connectome and NATMI. In matched cells it is closely related to the custom score by construction; these algorithms and their overlapping resource do not provide independent biological validation. The cell-type logFC score measures cell-type-versus-rest enrichment within a library, not disease gene log2 fold change. Magnitude consensus is transformed as 1-rank, so larger values indicate higher relative priority; its background candidate space may vary across animal libraries.

C, All four original illustrative pairs plus Col4a1–Itga3/Itgb1 are retained in a sign-and-eligibility audit. Each entry gives the sign of the MCAO-minus-Sham effect in discovery followed by validation. Custom entries require the cohort 10% expression gate. LIANA entries use 1-magnitude_rank and require complete observed animal scores. NA indicates unavailable eligibility or a nonfinite effect, including incomplete animal observation for LIANA; absence is not imputed as zero and does not prove absent signaling. A→E and E→A indicate the inferred astrocyte/endothelial directions. Colors redundantly encode the printed signs or NA status. No statistical significance is indicated by the color or sign.

Discovery contains 3 sham and 3 MCAO animal libraries; validation contains 3 sham and 2 MCAO libraries. Candidate denominators are not independent animal replicates. Exact zero effects have no positive/negative direction and are excluded from same-sign denominators; full finite-pair denominators are retained in source tables. No confidence interval treating edges as independent is shown. This figure reports methodological robustness and does not establish functional communication or barrier causation.
'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, default=BASE / 'tables')
    parser.add_argument('--output-base', type=Path, default=BASE)
    parser.add_argument('--export', action='store_true', help='Export final formats after inspecting the PNG preview')
    args = parser.parse_args()
    main(args.tables.resolve(), args.output_base.resolve(), args.export)
