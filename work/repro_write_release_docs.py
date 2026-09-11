"""Write release navigation and the twelve-point response from final audits."""
from pathlib import Path
import json,shutil
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2';OUT.mkdir(parents=True,exist_ok=True)
author=json.loads((ROOT/'repro_author.json').read_text(encoding='utf-8'))
record=OUT/'archive_record.json'
public=json.loads(record.read_text(encoding='utf-8')) if record.exists() else None
archive=('Published record: https://doi.org/'+public['doi']) if public and public.get('public_verified') else 'Public deposition is pending authenticated Zenodo access. Local archive metadata is not a DOI or proof of deposition.'
readme=f'''# Cross-cohort ligand–receptor reproducibility assessment

Current article: **Expression-level concordance and implementation sensitivity of ligand–receptor priority rankings in mouse cerebral ischemia**.

Author: {author['preferred_romanization']} ({author['name']}), {author['affiliation']}; ORCID {author['orcid']}. Correspondence is provided on the manuscript author page.

{archive}

## Read first

`manuscript_en.docx` and `manuscript_zh.docx` are English and Chinese versions of the same article. The English document is the international-submission working version; no journal submission is implied. The corresponding Markdown files preserve searchable text and linked references. `revision_response_12_points_zh.md` maps each requested revision to actual outputs and remaining publication requirements.

Five main figures cover samples/markers, the joint label baseline, aggregation implementation/ranking background, ambient-RNA sensitivity, and the added third cohort. Five supplementary figures cover identity, doublet scores, barrier-related transcription, paired sorted effects, and one spatial-expression map. Vector figures and their numerical sources accompany the PNG previews.

## Scientific scope and result families

- `null_*`: complete animal-label allocations, candidate keys, animal scores, observed effects, null distributions and all three original multiplicity families. The dynamic expression gate is recomputed per allocation. The later conditional 20/10 diagnostic remains separate.
- `rank_*`: native LIANA 1.10.0, the uniquely ranked score-column diagnostic, fixed-target percentiles, fixed full-network ranking universes, saturation/zero-effect and score-dependence audits. No installed package was patched. The known upstream fix was merged in [LIANA PR #261](https://github.com/scverse/liana/pull/261) on 2026-09-09. Its fixed-commit aggregation function exactly reproduced the unique-column diagnostic on all 22 raw networks; only this function, not a full newer release, was compared (`rank_upstream_*`). This case quantifies the correction's consequences and does not claim discovery of the defect.
- `ambient_*`: selected model-based corrected counts, marker summaries, >0 and >=1 custom eligibility sensitivity, and matched raw/corrected comparisons. Fractional estimated counts are not observed UMI counts and are never supplied to DESeq2 here.
- `third_*`: the GSE332910 raw primary-rule extension, author-reported three biological replicates/group, complete feature coverage, submitted library identities and matched score comparisons. Mouse or pool size per library remains unspecified. These descriptive results do not enter the original 200-allocation test.
- `sorted_*`: all 49,914 paired pool-level descriptive gene effects and pointwise intervals. Historical whole-transcriptome P/BH tables are preserved rather than discarded.

The main matched comparisons are on the same directed candidate sets. All-candidate and nonzero-effect concordance denominators are both retained. The joint null concerns labels in both original cohorts; it does not turn ten validation allocations into 200 validation animals. Native raw P values are not reused for diagnostic aggregation, corrected matrices or the third cohort.

## Archive layout

The code/results archive root contains `work/`, `outputs/`, `CITATION.cff`, `.zenodo.json`, `LICENSE` and `MANIFEST_SHA256.tsv`. The current release is `outputs/reproducibility_v2/`. Earlier tests and results remain labeled under `work/results`, `work/revision_analysis`, `work/sorted_rna` and earlier output directories; old manuscript drafts are not the entry point.

The code/results package contains integer-UMI pseudobulk counts, per-cell-type sample counts, estimated corrected aggregates, all statistical tests, full per-library ligand–receptor score networks and source/parameter/iteration audits. Large corrected single-cell NPZ matrices with cell/gene order are distributed separately. Original submitted GEO matrices remain downloadable using the archived accession manifests, URLs and checksums; they are not duplicated in the code archive. Pretrained reference models and original resources retain their source licenses and attribution.

## Reproduce the saved analysis

Use Python 3.12.14 and `work/requirements_reproducibility_lock.txt`. The effective analysis versions include NumPy 2.3.5, pandas 2.3.3, SciPy 1.18.1 and LIANA 1.10.0. Full package and actual-import records are in `work/repro_environment.json`. R fitting uses R 4.6.1, Bioconductor 3.23 and decontX 1.10.0, with complete per-run `sessionInfo.txt` and package inventories.

From a relocated release root:

```text
python work/repro_validate.py --require-complete
python work/repro_run_cached_pipeline.py
python work/repro_run_cached_pipeline.py --execute --include-completed-optional
```

The pipeline displays its plan unless `--execute` is supplied. It reconstructs tables from saved animal-level scores/aggregates, retains before-output copies and compares scientific content after reconstruction. It does not silently fit cell matrices or build the superseded manuscript. Optional `--figures` produces PNG previews that require visual review before vector export. See `work/repro_cached_pipeline_notes.md` for exact commands and audit meanings.

For a full raw-cell reconstruction, use the original/revision download and QC scripts, reference/doublet pipeline, and `work/repro_null_plan.md`, `repro_liana_rank_diagnostic_plan.md`, `repro_decontx_plan.md` and `repro_third_liana_plan.md`. Do not run all cell-matrix models concurrently on an 8 GB machine. The cached path is the practical first audit, while source scripts and count matrices preserve the full analytical chain.

## Numerical and source provenance

Default DecontX fits use maxIter=500. All and only fits failing the original convergence criterion receive a separately retained maxIter=2000 follow-up with unchanged input, labels, seed and model parameters. `work/repro_decontx_source.py` selects sources using this numerical rule, not communication results. Default and extended outputs, unresolved convergence flags, actual executed script snapshots and any audit-recovery records remain traceable.

The original 22 raw full-network diagnostic audits correspond to `work/repro_liana_rank_diagnostic_raw_run.py`. The current script adds corrected-source support and has a different hash. Original executed DecontX Python/R snapshots are retained. Lossless serialization changes are distinct from model changes. Where exact floating-point reconstruction matters for tied ranks, cached diagnostics use round-trip TSV parsing; integer and continuous estimates retain separate provenance.

## Publication status

Author identity, institution and correspondence have been supplied. Public archive verification and journal-specific submission declarations are separate from scientific execution. Funding, competing-interest declarations, author approval and any journal-specific AI-use statement must reflect the author's actual circumstances; this package does not invent them. A deposited code archive is not evidence that a manuscript has been submitted, reviewed or accepted.
'''
(OUT/'README.md').write_text(readme,encoding='utf-8')
response='''# 十二项修改落实情况

本文主线已从卒中候选机制挖掘改为有明确范围的计算可复现性研究。新增结果要求进一步修正原拟题：原生低排序一致性有一部分来自所用LIANA 1.10.0的可定位实现行为，因此最终标题讨论表达一致性和实现敏感性，不把原版本结果推广为所有配体—受体排序方法的普遍性质。

1. **归档。**代码、实际环境版本、细胞及pseudobulk汇总、全部旧检验、完整置换、网络中间分数、模型审计和图源表均纳入可搬移发布包；另保留大校正矩阵。作者为Sihuan Zhu（朱四欢）、安徽中医药大学，ORCID 0009-0005-9705-0585。真实公开记录状态见下文；未制造DOI。正文已移除用户指出的否定性仓库表述。
2. **定位。**标题、摘要、结果和讨论按可复现性问题重组。相同187条目标候选中，表达同向150/187（80.2%，ρ=0.725），原生优先级105/187（56.1%，ρ=0.228）。独特分数列仅排序一次后ρ=0.519，固定全局排名母集后ρ=0.645；同时完整报告零效应增加。改题由实际实现诊断决定。[上游PR #261](https://github.com/scverse/liana/pull/261)已于2026年9月9日合并同一修正，固定提交聚合函数在22份原始网络逐值吻合；本文贡献是其后果和剩余敏感性的定量案例，不声称首次发现该缺陷。
3. **Null baseline。**20×10=200种联合动物标签配置保持完整候选向量及共享配体/受体结构，动态10%门槛在每次重算。原79.9%对应零分布中位50.0%、95%参照24.5%–75.5%，精确P=0.005；6项探索统计Holm=0.030。固定集合30项及差值24项另行校正，最低分别0.15及0.12，全部结果保留。
4. **最小P。**增加了条件置换诊断：固定另一队列观测标签，发现/验证的动态同向率及ρ分别P=0.05/0.10。联合分辨率提高没有自动消除验证n=2，也不等于每队列分别证实复现；正文已说明这一差别。
5. **方法相关性。**表达与相对排序分两族。Connectome和NATMI共用expr_prod；expr_prod与lrscore完整网络名次近乎单调冗余。方法接口数不等于独立验证次数。诊断变体也不称官方修复版或新独立算法。
6. **环境RNA。**实际运行11个完整文库的DecontX，采用固定全脑参考分群、background=NULL。模型与数值收敛分别审计；默认未收敛者按统一上限2000追加，保留所有原结果。下游沿用固定细胞选择，以校正后全部基因总量归一化，重算表达及完整LIANA网络，并比较raw/corrected共同完整母集。>0与估计计数≥1的资格差别单列，校正失去候选不当作其原表达必属污染的证明。
7. **分选数据。**全部49,914项改为描述效应及点式95% t区间，图S4展示11个既定组分的全部6对池差值；不以全转录组显著数判断失败复现。与旧效应/区间数值误差<2×10⁻¹⁵，历史P/BH表原样保留。
8. **空间数据。**仅保留补充图S5的真实坐标表达；5个切片记录不称5只独立动物，空间相关不再承担主论证，其全描述表仍存档。
9. **更好的外部队列。**更新检索145条GDS命中并核实21个候选元数据，取得GSE332910及作者补充表S1，确认每组3个生物学重复、3个文库。6库实际下载、质控并分析；单核纹状体与原细胞数据差异及未分配核比例明确保留。在同一99条候选上，原两队列表达同向80/99，第三与它们为53/99和58/99（ρ=0.551/0.576）；同时变化的取材、测量和平台因素不能单独归因。第三队列不替换或隐藏旧n=2结果。
10. **摘要。**按现象、同母集量化对比和方法学含义重写。具体实现干预成为解释来源；联合P的适用范围明确，不以显著/不显著差异声称算法家族优劣。
11. **讨论。**讨论主体解释共用列重复排序、排名母集、RRA饱和与分母变化，结合去污染和第三队列说明范围；限制集中陈述。已纳入现有差异通讯及跨样本基准文献，不声称首次开展该类研究。
12. **图1C/D。**加入Kdr后白线及Target lineage/Myeloid组标注，保持与其他面板分组线一致；保留逐动物图源数据及矢量版本。

上述完成状态以最终validation、模型与图表审计为准。公开归档、代码重跑、文档排版和正式投稿是分别核验的步骤。

'''+archive+'\n'
(OUT/'revision_response_12_points_zh.md').write_text(response,encoding='utf-8')
shutil.copy2(OUT/'revision_response_12_points_zh.md',ROOT/'repro_revision_response.md')
print('Release README and twelve-point response updated')
