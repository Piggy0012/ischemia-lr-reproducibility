# Shared-structure null baseline: prospective computational amendment

Written 2026-09-11 before running the new null-baseline script or inspecting its outputs. This is a revision-stage analysis plan, not a preregistered original study. Existing disease contrasts and manuscript results were already known when this amendment was specified.

## Data and exchangeability unit

Use the original primary and reference_singlet selections, separately, in GSE174574 (three Sham, three MCAO) and GSE245386 (three Sham, two MCAO). Preserve each animal's entire multivariate ligand–receptor score vector. Exhaust all 20 × 10 = 200 joint condition-label allocations, including the observed allocation; never permute candidate edges, genes, or individual cells. The null assumes condition labels are exchangeable within each study under no treatment-label association. This does not remove unmeasured batch confounding, selection bias, or the validation cohort's n=2 limitation.

## Fixed common candidate universe (main comparison)

For each cell selection, use the intersection of directed Astrocyte–Endothelial resource edges reported with finite values for lr_means, expr_prod, lrscore, and magnitude_rank in every one of the 11 animal libraries. Require finite reconstructed custom scores for every library, too. This fixed set is determined without treatment labels; no group-dependent custom expression gate is added. Missing LIANA results are never imputed as zero. Reconstruct custom scores as the square root of sender minimum-complex mean-log-expression × receiver minimum-complex mean-log-expression from stored all-gene-normalized cell summaries. No new cellular computation is required for relabeling because normalization, cell selection, and LIANA within-library analysis did not use these disease labels.

Disease effect is MCAO mean minus Sham mean. Five score metrics are custom_score, lr_means, expr_prod, lrscore (expression-oriented family), and magnitude_priority = 1 − magnitude_rank (relative-ranking family). These are related computational summaries, not five independent validations. For each metric compare the two studies' complete vectors of disease effects using (1) same nonzero sign count divided by all fixed candidates, (2) the same count divided by candidates with nonzero effects in both studies, and (3) Spearman rho over all fixed candidates. Treat effects with absolute value ≤ 1e-12 as numerical zeros before computing signs and rho. Zero/zero is not counted as same direction. Report zero counts and both direction denominators.

For each of these three statistics, also compute each expression-oriented metric minus magnitude_priority. These gaps are evaluated under the same no-condition-label-association null. A tail probability for a gap is not a general test that the algorithms have equal reproducibility, and will not be described as such.

## Original group-gated custom-score sensitivity

Rebuild every resource edge having all constituent genes available in each cohort, irrespective of its original 5% screening status. For every joint label allocation independently, recompute the original 10% expression rule: at least two animals in either condition meet 10% minimum detection for all ligand subunits, and likewise for all receptor subunits. Intersect the two cohorts' eligible sets for that allocation, and calculate the same three custom-score concordance statistics. Report the varying candidate denominator. This repeats the selection step inside the null and provides a baseline for the original reported 323-edge comparison rather than conditioning on a disease-label-selected original set. A fixed original 323-edge null may be reported only as a clearly conditional sensitivity, not as the principal baseline.

## Tail probabilities and multiplicity

For each statistic, the one-sided exact upper-tail probability is the number of all 200 allocations having a statistic at least as large as observed (tolerance 1e-12), divided by 200. The observed allocation is included, so do not add a Monte Carlo +1 correction or generate artificial additional permutations. Report the null mean, median, and 2.5th/97.5th percentiles; these percentiles are a null reference interval, not an observed-effect confidence interval. Count finite allocations explicitly; abort if any intended statistic is undefined.

Three Holm adjustment families are fixed in advance: (A) 2 cell selections × 5 metrics × 3 concordance statistics = 30 tests; (B) 2 selections × 4 expression-minus-ranking gaps × 3 statistics = 24 tests; (C) 2 selections × 3 dynamically gated custom concordance statistics = 6 exploratory tests. Report every planned test, exact unadjusted P, and family-adjusted P; do not choose metrics, thresholds, or permutation tails based on the new results.

## Audit and follow-on interface

Save all 200 allocation statistics, observed summaries, full common-set per-animal scores and disease effects, fixed candidate keys, source-file hashes, allocation labels, script hash, and numerical checks. For a later mechanism decomposition only, export per-animal percentile ranks (average ties) computed over this same fixed candidate universe for each of the four expression-oriented metrics. These optional reranks are an interface, not an additional confirmatory null family in this amendment. Any later tests on them need a separate dated amendment.

Do not claim that a change from edge-level tests to an aggregate statistic automatically supplies power or significant results; conclusions depend on the observed exhaustive null distributions.
