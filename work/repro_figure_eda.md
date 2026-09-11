# Data profile: <DataFrame>

**Shape:** 30 rows × 17 cols

## Columns

| Column | Type | n | missing | summary |
|---|---|---|---|---|
| `config` | categorical | 30 | 0 | 2 levels: primary(15), reference_singlet(15); min_group_n=15 |
| `metric` | categorical | 30 | 0 | 5 levels: custom_score(6), lr_means(6), expr_prod(6), lrscore(6), magnitude_priority(6); min_group_n=6 |
| `statistic` | categorical | 30 | 0 | 3 levels: same_direction_all_fraction(10), same_direction_nonzero_fraction(10), spearman_rho(10); min_group_n=10 |
| `test_family` | categorical | 30 | 0 | 1 levels: A_fixed_common_metrics(30); min_group_n=30 |
| `observed` | continuous | 30 | 0 | mean=0.68, sd=0.159, range=[0.102, 0.807], skew=-2.28 (highly skewed); outliers=2 (IQR) |
| `null_mean` | continuous | 30 | 0 | mean=0.332, sd=0.239, range=[-4.44e-18, 0.5], skew=-0.67 (moderately skewed) |
| `null_median` | continuous | 30 | 0 | mean=0.332, sd=0.239, range=[0, 0.5], skew=-0.67 (moderately skewed) |
| `null_q025` | continuous | 30 | 0 | mean=-0.00703, sd=0.422, range=[-0.729, 0.373], skew=-0.76 (moderately skewed) |
| `null_q975` | continuous | 30 | 0 | mean=0.67, sd=0.122, range=[0.25, 0.776], skew=-2.14 (highly skewed); outliers=2 (IQR) |
| `n_allocations` | ordinal | 30 | 0 | 1 levels: 200(30); min_group_n=30 |
| `n_greater_or_equal_observed` | continuous | 30 | 0 | mean=9.9, sd=15.1, range=[1, 54], skew=1.74 (highly skewed); outliers=6 (IQR) |
| `exact_upper_tail_p` | continuous | 30 | 0 | mean=0.0495, sd=0.0755, range=[0.005, 0.27], skew=1.74 (highly skewed); outliers=6 (IQR) |
| `observed_n_candidates` | ordinal | 30 | 0 | 2 levels: 187(15), 174(15); min_group_n=15 |
| `null_n_candidates_min` | ordinal | 30 | 0 | 2 levels: 187(15), 174(15); min_group_n=15 |
| `null_n_candidates_max` | ordinal | 30 | 0 | 2 levels: 187(15), 174(15); min_group_n=15 |
| `holm_adjusted_p` | continuous | 30 | 0 | mean=0.348, sd=0.21, range=[0.15, 0.75], skew=1.22 (highly skewed); outliers=11 (IQR) |
| `n_tests_in_family` | ordinal | 30 | 0 | 1 levels: 30(30); min_group_n=30 |

## Group structure
- Grouped by: `config`, `statistic`
- Number of groups: 6
- Group size: min=5, median=5, max=5
- **WARN**: at least one group has n<10 — use box/violin + stripplot rather than mean-only bar chart.

## Correlations (Pearson, sorted by |r|)
- `null_mean` ↔ `null_median` : r = 1.000 (very strong)
- `n_greater_or_equal_observed` ↔ `exact_upper_tail_p` : r = 1.000 (very strong)
- `observed` ↔ `null_q975` : r = 0.974 (very strong)
- `null_mean` ↔ `null_q025` : r = 0.971 (very strong)
- `null_median` ↔ `null_q025` : r = 0.971 (very strong)
- `exact_upper_tail_p` ↔ `holm_adjusted_p` : r = 0.934 (very strong)
- `n_greater_or_equal_observed` ↔ `holm_adjusted_p` : r = 0.934 (very strong)
- `observed` ↔ `n_greater_or_equal_observed` : r = -0.789 (very strong)
- `observed` ↔ `exact_upper_tail_p` : r = -0.789 (very strong)
- `observed` ↔ `holm_adjusted_p` : r = -0.740 (very strong)
- ... +18 more pairs

## Warnings
- 列 'metric' 至少有一个类别 n<10 — 小样本必须展示原始数据点，不要只画均值柱状图。

## Chart suggestions (preliminary)
- 分类 vs 连续，小样本（每组 n<10）→ **箱线图/小提琴图 + stripplot 叠加原始点**；**避免**只画均值柱状图，会掩盖分布。
- ≥3 个连续变量 → 相关性热力图（['observed', 'null_mean', 'null_median', 'null_q025', 'null_q975']）或 pairplot 散点矩阵
- 分类维度组合数 = 240（config, metric, statistic, test_family, n_allocations, observed_n_candidates, null_n_candidates_min, null_n_candidates_max, n_tests_in_family 全交叉），**一张图塞不下**——建议按某一维拆成多面板，或选择子集。
- observed 高度偏态（skew=-2.28）→ 考虑对数变换或小提琴图代替均值柱图
- null_q975 高度偏态（skew=-2.14）→ 考虑对数变换或小提琴图代替均值柱图
- n_greater_or_equal_observed 高度偏态（skew=1.74）→ 考虑对数变换或小提琴图代替均值柱图
- exact_upper_tail_p 高度偏态（skew=1.74）→ 考虑对数变换或小提琴图代替均值柱图
- holm_adjusted_p 高度偏态（skew=1.22）→ 考虑对数变换或小提琴图代替均值柱图

> 这是基于数据形态的**初步建议**。最终图型选择必须结合**论证目标**（你想说什么）—— 详见 `references/chart_selection.md`。