# Exploratory third-cohort full-network implementation sensitivity

Recorded on 2026-09-11 before running LIANA on GSE332910. The first custom
coavailability comparison is already known: 63.23% / rho 0.560 against
GSE174574, and 60.00% / rho 0.547 against GSE245386, using pairwise 10%-eligible
candidate sets. This follow-on analysis is therefore exploratory.

Use the six submitted RNA libraries with explicit Sham or 24-hour titles.
Supplementary Table S1 reports three biological replicates and three libraries
per condition. The individual mouse or pool membership is not supplied. Keep
library-level identifiers and report author-reported biological replication;
do not silently assert one mouse per library or add these libraries to the
original two-cohort label-permutation analysis. The comparison is descriptive,
without additional P values.

Retain the original primary marker labels and QC thresholds, without tuning
to third-cohort effects. In each library, retain assigned classes with at least
30 cells. Stream the gene-major matrix into a context-only sparse matrix,
keeping every submitted gene and the original full-cell UMI denominators.
Verify the retained cells' UMI/detected-gene totals against QC and each retained
class's pseudobulk against the independent stream aggregation.

Apply LIANA 1.10.0 with the same five-method default, mouseconsensus resource,
float32 log1p(CP10k), expr_prop 0.1, min_cells 30, n_perms 100, seed 20260910,
n_jobs 1, use_raw False, return_all_lrs False and full cell-type network as the
original two-cohort full-network diagnostic. Preserve all native full-network
scores and eligible target rows. The third cohort has no prior LIANA export;
do not create a false comparison to old scores. Instead, reconstruct its native
magnitude rank from the full-network scores and require agreement within 1e-12.
Then rank each unique magnitude score column once on exactly the same network
using the installed rank aggregation function. Label this output as an
implementation diagnostic, not an official package correction or independent
biological validation. Record exact score specifications and source hashes.

Compare six metrics: custom coavailability, lr_means, expr_prod, lrscore,
native magnitude priority (1 minus rank), and unique-column diagnostic priority.
The first four are expression-related summaries, not four independent tests.
Estimate each disease effect as the equal-library MCAO mean minus Sham mean.
Report direction agreement with and without zero effects, zero counts, Spearman
rho and exact rank ties. Effects within 1e-12 of zero are treated as zero.

For each older cohort versus the third, use an identical candidate set across
all six metrics, requiring finite native/diagnostic/method scores in all
libraries of the pair (12 and 11 libraries, respectively), with no imputation of
ineligible edges. Also export both comparisons restricted to the same candidate
set complete in all 17 libraries across the three cohorts. This distinguishes
pairwise transfer from changes caused by candidate coverage. These universes
were specified before the third-cohort LIANA results. Full-network ranking is
computed before restricting to any target or cross-cohort candidate subset.
Persisted numerical TSVs are parsed with pandas `float_precision='round_trip'`
to preserve their serialized double-precision values and exact ties during
summary and any subsequent percentile/rank calculations.

Build only one sparse context cache at a time. Cache building is measured below
475 MB peak working set. Do not run LIANA concurrently with the larger DecontX
jobs; coordinate resource availability first. Stop before another library if
the LIANA process's observed peak working set exceeds 2 GiB.
