# Cached v2 validation and reconstruction

The release root contains `work/` and `outputs/`. Both entry points resolve paths from their own `__file__`, so the directory can be relocated without preserving the original Windows username or checkout. Do not use an earlier manuscript pipeline as the v2 entry point.

## Environment

Use the archived `requirements_reproducibility_lock.txt`, `requirements_revision.txt`, and software-version records. The verified analysis environment uses Python with NumPy 2.3.5, pandas 2.3.3, SciPy 1.18.1, LIANA 1.10.0, statsmodels, matplotlib, Pillow and the transitive LIANA dependencies. The raw-run script snapshot is `repro_liana_rank_diagnostic_raw_run.py`; its hash is compared against all 22 original full-network audits. The current diagnostic script additionally supports corrected-count source selection and is a different, intentionally retained file.

Validation reads small tabular outputs, compressed summaries and JSON audits. It verifies input hashes, frozen candidate counts, the complete animal score matrices, equal-animal disease effects, exact allocation tails and independently recomputed Holm corrections within the original families. It retains the distinction between original animal-resolved libraries, third-cohort author-reported biological replicates with unspecified mice per library, and six matched sorted-sample pool comparisons. Unknown animal counts are missing metadata, not zero animals. Native within-library LIANA P values and descriptive diagnostics are not promoted to disease-level or adjusted inference.

## Commands

From a relocated release root, using the installed analysis Python:

```text
python work/repro_validate.py
python work/repro_validate.py --require-complete
python work/repro_run_cached_pipeline.py
python work/repro_run_cached_pipeline.py --execute
python work/repro_run_cached_pipeline.py --execute --include-completed-optional
python work/repro_run_cached_pipeline.py --execute --figures
```

The first validation command checks completed stages and explicitly lists pending stages. `--require-complete` exits nonzero until ambient-RNA and third-cohort LIANA completion evidence exists and passes. Exit code 1 indicates failed scientific/provenance checks; exit code 2 indicates otherwise valid results with required stages still pending. Reports are written to `outputs/reproducibility_v2/validation/` by default. Validation never modifies source tables, original tests or model outputs.

The cached pipeline prints its plan unless `--execute` is supplied. Default stages rebuild null, fixed-target rank, raw target summary, the isolated fixed-commit upstream correction audit, fixed-global background, conditional-null, sorted descriptive and third custom-expression tables. The optional switch only includes already completed ambient/third-LIANA summaries; absent stages are explicitly skipped and cannot trigger model fitting. The `--stages` option permits a bounded subset in the order supplied; selecting a downstream stage does not silently recreate missing prerequisites. A full default run provides the intended dependency order.

Each execution creates a new `work/repro_cached_runs/<UTC-run-id>/` directory containing the plan, full before-output copies, before/after validation, protected-input hashes, per-stage logs, a semantic comparison table and a run audit. No previous run directory is overwritten. Source historical tests, animal aggregates and full-network raw score exports are protected. Differences in table schemas, row order, candidate numbers, missingness, strings or numerical values beyond the reported tolerance fail the run. Gzip timestamps can change despite identical decoded tables; byte hashes and semantic comparison are therefore both reported. JSON provenance changes remain visible rather than being removed from the audit. Validation and comparison use pandas `float_precision='round_trip'`: the default parser can perturb nearly tied percentile effects and change a recalculated Spearman correlation despite adequate decimal serialization. No broader tolerance is substituted for exact restoration of these stored floats.

## Cache and model boundaries

The original two-cohort expression null reconstructs scores from saved per-animal mean-log-expression and detection-fraction summaries. It does not scan raw cells. Fixed-global ranking reads already saved full-network scores and calls only the aggregation function. Raw/corrected target-summary stages call `summarize()` directly; the inference entry points are disabled in the worker. Third custom reconstruction requires all six existing aggregation audits, means, fractions and library-count summaries; the streaming fallback is disabled. A worker also denies matrix-format reads and network connections, so missing caches fail visibly.

The large raw GEO and corrected cell-count matrices, DecontX fitting, CellTypist/Scrublet inference and full-network LIANA computation are separate upstream stages. They are never run by this entry point. Preserve their download checksums, model/source audits and supplied model-version records. Additional raw/corrected matrix archives may be required to reproduce those upstream fits, while the code/results archive supports the explicitly bounded cached reconstruction.

## Figures and incomplete model stages

By default, table reconstruction does not rewrite figures. `--figures` generates Figure 1/2/3 and sorted-data PNG previews where their selected stage supports them. It does not export publication PDF/SVG or claim human visual inspection. The old exports and original review evidence are preserved in the before-output directory; newly generated PNGs require actual visual inspection before their corresponding figure scripts are invoked with `--export`. Review includes matching source values, both concordance denominators, zero-effect counts, readable text and embedded PDF fonts. A numeric check cannot substitute for this review.

Ambient completion means that all planned libraries and numerical follow-up decisions have finished and are documented. The resolver selects the default 500-iteration fit if its recorded convergence threshold was reached, otherwise the separately retained uniform 2,000-iteration follow-up. Original fits remain immutable. Selected source audit hashes and the matched raw/corrected all-animal/all-metric candidate audit must accompany the final sensitivity outputs. A finished model fit is not automatically a converged fit; unresolved convergence is disclosed. No native raw null P value is transferred to corrected, unique-column, fixed-global or third-cohort diagnostics.

The third cohort reports three biological replicates and three RNA libraries per group in Supplementary Table S1. Single-mouse versus pool size remains unspecified. Its added analyses stay descriptive and outside the original 200 joint allocations. The conditional 20/10 label diagnostic also remains separate from the original A/B/C multiplicity families and does not remove the validation cohort's n=2 limitation.

## Completed relocation check, 2026-09-11

The relocated release at `work/repro_delivery/analysis` passed all 12 cached stages and all 2,749 complete validation checks, with no protected-input changes, scientific table changes or pending stages. The successful local run is `work/repro_cached_runs/20260911T065221_225850Z` relative to that relocated release. Compact evidence is retained in `outputs/reproducibility_v2/validation/relocated_cached_pipeline` in both the original and relocated directories. The same frozen interpreter was used; this verifies directory relocation and cached source resolution, not a fresh-machine dependency installation.

Two earlier failures are deliberately retained. `validation/relocated_initial` records an older gzip header inherited for the sorted normalized source despite identical decompressed content; staging now copies the exact current audited bytes. `validation/relocated_path_failure` records a complete numerical reconstruction followed by a validator path error when an outer `work/` directory preceded the relocated `outputs/` directory. The resolver now selects the innermost archive-root marker and still verifies exact source hashes. The complete rerun passed without changing scientific definitions or relaxing numeric checks. Detailed metadata differences identify only regenerated gzip hashes, relocated provenance paths and measured process-memory changes. Local duplicate before-output backups remain available but are excluded from the release ZIP in favor of the complete compact audit evidence.
