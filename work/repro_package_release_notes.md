# Release packaging

`repro_package_release.py` packages an existing finalized stage. Its default invocation is a read-only plan; it never runs the staging script, fits a model, downloads a dataset or publishes an archive.

```text
python work/repro_package_release.py --with-corrected-matrices
python work/repro_package_release.py --execute --with-corrected-matrices
```

The default stage is `work/repro_delivery/analysis`, produced separately by `repro_stage_archive.py`. Output is a sibling `work/repro_delivery/release` directory, so the ZIPs cannot include themselves. A nonempty release output directory is never overwritten; use a new explicitly named directory for a subsequent version. Source locations are resolved relative to `__file__`, and stage/output locations must remain in the declared workspace.

Execution first runs the current staged `repro_validate.py --require-complete` with the analysis interpreter. The validator source must match the current reviewed copy. Any failed check or pending ambient/third-LIANA stage prevents packaging. The new complete audit and log are included in the staged code/results payload. A partial model status cannot become a complete release by renaming a file or reusing an earlier partial validation report.

`code_results.zip` contains the staged analysis code, environment records, derived counts/scores, all original and revision tests, results, figures and source audits. It excludes interpreter caches, self-referencing manifest/integrity files, and `work/repro_cached_runs` (local duplicate before-output backups). Successful same-workspace and relocated cached reconstruction evidence is preserved under `outputs/reproducibility_v2/validation/`, including plans, logs, before/after validations, protected-source hashes and semantic comparisons. Local backup directories remain available and are never deleted by packaging. Runtime directories, account files, keys, external symlinks and nested archives are rejected. Sources are restricted to this explicit stage; the script does not inspect or collect unrelated computer files.

The optional `corrected_matrices.zip` contains all 11 corrected count matrices selected by `repro_decontx_source.resolve`, each with its matching gene and cell order and completed audit. For a default 500-iteration fit selected for the uniform 2,000-iteration numerical extension, the original default matrix and metadata are also retained at their original, distinct path. Archive paths preserve `work/repro_decontx/...` and `work/repro_decontx_maxiter2000/...`, allowing extraction beside the code/results archive. NPZ is retained and the duplicate Matrix Market representation is omitted. Matrix hashes must match their per-library model audit, and selected/default audit hashes must match the code stage's `ambient_selected_sources.tsv`. The number of selected libraries remains 11 even when additional original nonconvergent matrices are archived.

Each ZIP includes `MANIFEST_SHA256.tsv`, listing every payload member, byte count, SHA256, role and compression policy. An identical external `<archive>_SHA256.tsv` is supplied. The manifest does not contain its own checksum. External `release_integrity.json` records the ZIP and manifest checksums, full member-count and CRC checks, member SHA256 verification, complete-validation audit checksum and convergence source decisions; it is not inserted into an archive that it hashes. This avoids circular hash dependencies.

Files already compressed as gzip, NPZ, PNG, PDF, Office documents and R data are stored directly; other files use DEFLATE level 1. Every member is checked with ZIP CRC and with the full SHA256 manifest. Source hashes are rechecked after compression to detect concurrent changes. Final ZIP names are assigned only after these checks pass; a failed partial archive remains clearly suffixed `.partial` for diagnosis.

Packaging verifies integrity and provenance. It does not create a DOI or claim public deposition. The final repository/Zenodo identifier must be recorded after the separate publication step.

## Staging review as of 2026-09-11

The current staging script's `repro_*.py/md/json/R` patterns cover the validator, guarded cached pipeline, package entry point, English source drafts, numerical-selection resolver and raw-run script snapshot. Its explicit directory list now includes raw/corrected full LIANA networks, the new third-cohort LIANA outputs, third-cohort mean/fraction/count caches and the third context-cache metadata. The default and maxIter=2000 model metadata are both included, while corrected matrices remain in the separate package. Existing original derived sources are inherited from the audited earlier stage.

The stage was subsequently built and its initial validation exposed one stale compressed source inherited from the earlier release. The normalized sorted-expression values were identical after decompression, but its gzip SHA256 differed from the current audit. Staging now copies the exact current audited source; the original failure, both hashes and the content-equivalence check are retained in `validation/relocated_initial/`. Complete staged validation passed 2,749 checks after this correction. Exact historical Python and R source snapshots and their provenance maps are retained separately from current wrappers.

The complete 12-stage relocated reconstruction then exposed and resolved an inner/outer archive path-marker ambiguity in the validator. The corrected validator passed the final 2,749 checks with no scientific or protected-source changes. The failed-path evidence and final successful evidence are retained under `validation/relocated_path_failure` and `validation/relocated_cached_pipeline`. The packaging dry run found all 11 convergence-selected model sources complete. Publication packaging must still run its fresh complete validation after final manuscript and visual-QA synchronization.
