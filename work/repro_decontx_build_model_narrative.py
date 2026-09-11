"""Summarize frozen final DecontX sources using small tables and NPZ headers only."""
from pathlib import Path
import hashlib
import json
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    default = json.loads((ROOT / 'repro_decontx/completion_status.json').read_text())
    extension = json.loads((ROOT / 'repro_decontx_maxiter2000/completion.json').read_text())
    assert default['completed_samples'] == default['expected_samples'] == 11
    assert extension['all_completed'] and extension['all_refits_converged']
    source_integrity = json.loads((ROOT / 'repro_decontx/source_integrity_check.json').read_text())
    assert source_integrity['all_declared_source_hashes_match']
    assert source_integrity['completed_outputs_checked'] == 12
    records, tables = [], []
    for sample in default['samples']:
        selected = 'repro_decontx' if sample['convergence_threshold_reached'] else 'repro_decontx_maxiter2000'
        folder = ROOT / selected / sample['dataset'] / sample['sample']
        audit_path = folder / 'audit.json'
        a = json.loads(audit_path.read_text())
        assert a['convergence_threshold_reached']
        c = pd.read_csv(folder / 'contamination.tsv.gz', sep='\t', usecols=['contamination', 'raw_total', 'corrected_total'])
        assert len(c) == a['n_cells']
        assert np.isfinite(c.to_numpy()).all() and (c.to_numpy() >= 0).all()
        assert c.contamination.le(1).all()
        with zipfile.ZipFile(folder / 'corrected_counts.npz') as archive:
            with archive.open('data.npy') as handle:
                version = np.lib.format.read_magic(handle)
                assert version == (1, 0)
                shape, fortran_order, dtype = np.lib.format.read_array_header_1_0(handle)
        assert dtype == np.dtype('float64')
        assert int(c.corrected_total.eq(0).sum()) == a['n_zero_corrected_total']
        assert np.isclose(c.raw_total.sum(), a['raw_total'], rtol=0, atol=1e-7)
        assert np.isclose(c.corrected_total.sum(), a['corrected_total'], rtol=1e-12, atol=1e-5)
        records.append({
            'dataset': sample['dataset'], 'sample': sample['sample'],
            'condition': sample['condition'], 'selected_source': selected,
            'audit_relative_path': str(audit_path.relative_to(ROOT)).replace('\\', '/'),
            'audit_sha256': sha(audit_path),
            'contamination_sha256': sha(folder / 'contamination.tsv.gz'),
            'n_cells': len(c), 'maxIter': a['requested_parameters']['maxIter'],
            'last_iteration': a['last_logged_iteration'],
            'last_max_change': a['last_logged_max_divergence'],
            'converged': a['convergence_threshold_reached'],
            'median_estimated_contamination': float(c.contamination.median()),
            'mean_estimated_contamination': float(c.contamination.mean()),
            'raw_total': int(c.raw_total.sum()),
            'corrected_total': float(c.corrected_total.sum()),
            'n_zero_corrected_total': int(c.corrected_total.eq(0).sum()),
            'data_dtype': str(dtype), 'fractional_values': a['fractional_values'],
        })
        tables.append(c)
    pooled = pd.concat(tables, ignore_index=True)
    report = {
        'n_libraries': len(records), 'n_cells': len(pooled),
        'default_converged': sum(s['convergence_threshold_reached'] for s in default['samples']),
        'final_selected_all_converged': all(r['converged'] for r in records),
        'pooled_cell_median_estimated_contamination': float(pooled.contamination.median()),
        'pooled_cell_mean_estimated_contamination': float(pooled.contamination.mean()),
        'library_median_estimated_contamination_min': min(r['median_estimated_contamination'] for r in records),
        'library_median_estimated_contamination_max': max(r['median_estimated_contamination'] for r in records),
        'n_zero_corrected_total': int(pooled.corrected_total.eq(0).sum()),
        'raw_total': int(pooled.raw_total.sum()),
        'corrected_total': float(pooled.corrected_total.sum()),
        'removed_total_fraction': float(1 - pooled.corrected_total.sum() / pooled.raw_total.sum()),
        'count_data_dtype': 'float64',
        'pooled_summary_is_descriptive_not_replication': True,
        'extended_sample_comparison': extension['comparisons'],
    }
    assert report['n_cells'] == 115355 and report['default_converged'] == 10
    assert report['n_zero_corrected_total'] == 0
    pct = 100 * report['pooled_cell_median_estimated_contamination']
    low = 100 * report['library_median_estimated_contamination_min']
    high = 100 * report['library_median_estimated_contamination_max']
    methods_zh = (
        '将 DecontX 作为计算性环境 RNA 敏感性分析逐库运行（R 4.6.1、Bioconductor 3.23、decontX 1.10.0）。'
        '输入为 GSE174574 和 GSE245386 的 11 个文库中全部原始质控合格细胞的整数 UMI 计数；重复基因符号在模型拟合前合并求和。'
        '固定使用独立全脑参考注释的 whole_brain_broad 标签作为 z，不按目标候选或拟合结果重新选择细胞。'
        '参数为 background=NULL、batch=NULL、seed=20260911、delta=c(10,10)、estimateDelta=TRUE、convergence=0.001、iterLogLik=10，数值线程设为 1。'
        '初始 maxIter=500；保留全部初始输出，并按统一数值规则仅对未达到收敛阈值的文库用相同输入、分群和随机种子将 maxIter 延长至 2000 重新拟合。'
        '最终来源仅依据数值收敛选择：初始收敛者使用默认输出，未收敛者使用延长拟合输出；选择不依据候选通讯结果。'
        '校正矩阵保存为 float64 的分数估计，提供 cells×genes CSR NPZ 与 genes×cells MatrixMarket；不取整、不将其作为原始 UMI 或输入 DESeq2，并保留基因及细胞顺序、校正后每细胞总量和来源 SHA256。'
        '该模型在缺少空液滴实测背景时从已过滤细胞推断非本细胞来源成分，依赖所给细胞群标签；参考身份不确定及细胞群间共享的真实表达可能影响分离结果。'
        '污染比例因此属于模型估计而非实测污染率；该分析不直接测量配体释放、受体活化或屏障功能。'
    )
    methods_en = (
        'DecontX was run separately for each library as a computational ambient-RNA sensitivity analysis (R 4.6.1, Bioconductor 3.23, decontX 1.10.0). '
        'Inputs were integer UMI counts for all originally QC-passing cells in 11 libraries from GSE174574 and GSE245386; counts for duplicate gene symbols were summed before fitting. '
        'The fixed whole_brain_broad labels from an independent whole-brain reference annotation were supplied as z, without reselecting cells according to target candidates or fitted results. '
        'Parameters were background=NULL, batch=NULL, seed=20260911, delta=c(10,10), estimateDelta=TRUE, convergence=0.001, and iterLogLik=10, with numerical threads set to one. '
        'The initial iteration limit was 500. All initial outputs were retained, and a uniform numerical rule required only libraries failing the convergence threshold to be refitted with maxIter=2000 using identical inputs, groups, and seed. '
        'Final sources were selected solely by this convergence rule: default outputs for initially converged libraries and extended-fit outputs for initially unconverged libraries, irrespective of candidate communication results. '
        'Corrected values were retained as fractional float64 estimates in cells-by-genes CSR NPZ and genes-by-cells MatrixMarket formats, with preserved cell and gene order, corrected per-cell totals, and SHA256 provenance. '
        'They were neither rounded nor treated as observed UMI counts or supplied to DESeq2. '
        'Without a measured empty-droplet background, the model infers non-native components from filtered cells conditional on the supplied groups; uncertain reference identities and genuinely shared expression across groups may affect this separation. '
        'Contamination fractions are therefore model estimates rather than measured contamination rates, and this analysis does not directly measure ligand secretion, receptor activation, or barrier function.'
    )
    results_zh = (
        f'完成了 11 个文库共 115,355 个原始质控合格细胞的 DecontX 拟合。10 个文库在默认 500 次迭代上限内达到 0.001 收敛阈值；唯一未收敛文库 GSM5319992 在相同输入和设置、仅将上限改为 2000 后于第 637 次迭代达到阈值（末次最大参数变化 0.0009979）。'
        f'最终采用 10 份默认结果及该文库的延长拟合结果，11 个最终来源均达到预设数值阈值，全部默认结果仍保留。'
        f'合并细胞后的估计污染比例中位数为 {pct:.2f}%，各文库中位数范围为 {low:.2f}%–{high:.2f}%；合并统计仅作描述，不把细胞作为独立生物学重复。'
        '校正后每细胞总量均大于零，11 个最终矩阵均为 float64 分数估计。'
        'GSM5319992 延长前后估计污染中位数由 6.32% 变为 6.14%，计数矩阵绝对差总和相当于其原始总 UMI 的 0.3565%；不过单细胞污染估计最大绝对差为 37.54 个百分点，故总体差异较小不能替代逐细胞或下游结果核查。'
        '这些结果描述模型拟合与数值敏感性，不单独支持特定通讯候选或屏障调控机制。'
    )
    results_en = (
        f'DecontX fitting was completed for 115,355 originally QC-passing cells in 11 libraries. Ten libraries reached the 0.001 convergence threshold within the default 500-iteration limit. '
        f'The only initially unconverged library, GSM5319992, reached the threshold at iteration 637 after refitting identical inputs and settings with only the iteration limit increased to 2000 (last maximum parameter change, 0.0009979). '
        f'The final sources comprised ten default fits and this extended fit; all eleven met the specified numerical threshold, and every default result was retained. '
        f'The pooled cell-level median estimated contamination fraction was {pct:.2f}%, and library-specific medians ranged from {low:.2f}% to {high:.2f}%. '
        'Pooled summaries were descriptive and did not treat cells as independent biological replicates. '
        'All cells retained positive corrected totals, and all eleven final matrices contained fractional float64 estimates. '
        'For GSM5319992, extending the fit changed the median estimated contamination from 6.32% to 6.14%; the sum of absolute matrix differences was 0.3565% of its original total UMI count. '
        'However, the maximum absolute change in an individual cell\'s estimated contamination was 37.54 percentage points, so a small aggregate change cannot replace cell-level or downstream checks. '
        'These results characterize model fitting and numerical sensitivity and do not by themselves support particular communication candidates or barrier-regulatory mechanisms.'
    )
    result = {
        'status': 'complete', 'methods_zh': methods_zh, 'methods_en': methods_en,
        'results_zh': results_zh, 'results_en': results_en,
        'statistics': report, 'final_source_selection': records,
        'references': ['https://bioconductor.org/packages/release/bioc/html/decontX.html',
                       'https://doi.org/10.1186/s13059-020-1950-6'],
        'summary_script_sha256': sha(Path(__file__)),
        'source_audits_modified': False,
    }
    (ROOT / 'repro_decontx_model_narrative.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    pd.DataFrame(records).to_csv(ROOT / 'repro_decontx/final_selected_source_summary.tsv', sep='\t', index=False)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
