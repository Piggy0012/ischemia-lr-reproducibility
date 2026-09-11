# Upstream correction: attribution and scope

Verified on 2026-09-11 from the official pull-request API and fixed merge source. [PR #261](https://github.com/scverse/liana/pull/261) was merged at 2026-09-09T07:06:26Z as [d4211373692e7b9c10210488ccb1efe06452b097](https://github.com/scverse/liana/commit/d4211373692e7b9c10210488ccb1efe06452b097). The aggregation source is pinned to [the fixed merge, new `_core` path](https://github.com/scverse/liana/blob/d4211373692e7b9c10210488ccb1efe06452b097/src/liana/_core/_pipe_utils/_aggregate.py). Its SHA256 is `484b2e61d7bb01d141116a6fdf8f8b7d269973396be3354d10702222770aa875`. API responses, exact source, executable AST-identical function subset and license are retained in `work/repro_literature/`.

## 中文方法补充

在最终来源复核中，确认上游PR #261已于2026年9月9日合并针对共享评分列重复排名的修正。我们将固定合并提交中的4个聚合函数单独提取，逐个函数核对AST一致性，在现有NumPy/pandas/SciPy环境中执行；未安装新的LIANA包。对全部22个原始完整网络恢复3个表达评分列的float32精度，保留原网络、参数与187/174条目标候选，仅重算该固定提交的聚合函数。另分别以相同float64输入复核两个聚合实现，记录排名精度而不改变表达顺序。

## 中文结果补充

固定上游聚合函数与此前“独特评分列各排名一次”的诊断在全部22个完整网络及各自目标集合上逐值完全一致，最大绝对误差为0。在本环境中，SciPy 1.18.1的rankdata对float32输入保留float32，因此两实现的排名矩阵精度相同；两者接受相同float64输入时亦逐值一致。跨队列ρ仍为0.519/0.514；全部候选同向数仍为122/187和115/174，非零效应分母为168/158。该追加核验不产生新的P值，未更改原零分布或多重校正家族。

## 中文讨论归属修订

共享expr_prod列的重复排名问题及其修正应归于已合并的上游工作。本研究的增量是量化这一已知聚合修正在相同动物、评分和候选条件下对跨队列结果的影响，并继续评价候选母集与评分饱和带来的残余敏感性。官方PR还修改了表达比例等默认设置；本次核验仅执行固定提交的聚合函数，因此不能表述为完整新版LIANA流程的比较。稿件应删除“首次发现该bug”或暗示上游尚未修正的表述，同时保留明确的软件版本与函数级审计范围。

## English Methods insertion

Final source verification identified an already merged upstream correction in LIANA PR #261, merged on 9 September 2026. Four aggregation functions were extracted from the fixed merge commit d4211373692e7b9c10210488ccb1efe06452b097, with function bodies verified as AST-identical. These functions were executed in the existing numerical environment on all 22 saved raw complete-network score tables, restoring the three magnitude inputs to float32 and preserving network membership, method specifications and the 187/174 fixed target candidates. No newer LIANA package was installed. A separate precision check supplied both implementations with identical float64 copies of the same score values.

## English Results insertion

The fixed upstream aggregation function and our earlier unique-column diagnostic were bitwise identical across all 22 complete networks and their fixed target subsets, with a maximum absolute difference of zero. In the audited SciPy 1.18.1 environment, rankdata retained float32 precision for these inputs; the two implementations were also identical when both received float64 inputs. Cross-cohort correlations remained 0.519/0.514, with concordance counts of 122/187 and 115/174 overall and nonzero-effect denominators of 168/158. The additional audit introduced no P values and left the original null families unchanged.

## English Discussion replacement

The shared-column ranking issue and its correction were already documented upstream. Our contribution is to quantify the cross-cohort consequences of that known correction while preserving animal inputs, expression scores and candidate membership, and to examine the remaining sensitivity to the ranking universe and score saturation. The merged PR also changed expression-proportion defaults and other pipeline behavior. The present equivalence audit therefore applies to the isolated aggregation functions under the retained inputs and environment, rather than to a complete newer-version LIANA analysis. Attribution to the upstream correction replaces any implication that this study discovered a previously unknown software defect.

## Source-data and audit record

`rank_upstream_per_library_equality.tsv` contains both complete-network and target-subset checks. `rank_upstream_animal_scores.tsv.gz`, `rank_upstream_effects.tsv.gz`, and `rank_upstream_summary.tsv` retain the descriptive quantities. `rank_upstream_audit.json` records all input/source hashes, AST extraction, software versions, the absence of new inference, and observed memory below the 400 MiB bound. The initial isolated execution peaked at 177,618,944 bytes (169.4 MiB); the later guarded cached reconstruction peaked at 204,165,120 bytes (194.7 MiB). The original raw/native, unique-column and fixed-global outputs remain intact.
