# Lossless native score serialization

2026-09-11, after the first complete native model and before any cross-cohort CellChat effect comparison.

The first library's exact native RDS and decimal `write.table` output differed in 2,684 probability values, by at most approximately 5e-16 absolute. R `write.table` does not necessarily preserve 17 significant digits merely because `options(digits=17)` is set. The number of distinct probability values was unchanged in this probe, but downstream ranks must not rely on presumed harmless rounding.

No native model, seed, cohort label, eligibility rule or candidate set changes. Preserve the original full-precision RDS and original TSV, including the originally recorded model/output hashes. Export an additional `full_network_lossless.tsv.gz` directly from each native RDS using `sprintf('%.17g')` before writing numeric fields as characters. Require exact equality of every probability and internal P value after an independent R decimal readback. Record hashes of source RDS, original TSV, new TSV and the exporter. Python summaries read this new file using pandas `float_precision='round_trip'`.

This is a persistence correction with a complete provenance chain. All source × target × ligand/receptor rows, zero scores and P values remain present. The running native inference script is not edited or rerun to make this correction.
