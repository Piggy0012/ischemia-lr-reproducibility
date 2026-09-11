# Actual DecontX sensitivity analysis

The amendment `../repro_decontx_plan.md` was written before outputs were examined. Only libraries with a final `audit.json` are complete. Input raw count files remain unchanged.

The default commands run from the study workspace (Windows PowerShell):

```powershell
& 'work/r_runtime/bin/Rscript.exe' --vanilla 'work/repro_decontx_install.R' 'work/repro_decontx/rlibrary'
& 'work/revision_env/Scripts/python.exe' 'work/repro_decontx.py' GSE174574 GSE245386
```

R 4.6.1 was downloaded from the official CRAN URL and its MD5 verified. Bioconductor 3.23/decontX 1.10.0 was installed from official binary repositories. See `installers/r_runtime_manifest.json`, package installation logs, and each sample's `sessionInfo.txt` / `installed_packages.tsv`. Runtime and dependency folders are local execution infrastructure, not necessary study-result deliverables.

For another machine, provide an appropriate Rscript and package library using `DECONTX_RSCRIPT` and `DECONTX_RLIBRARY`. To regenerate into a separate location without altering audited results, set `DECONTX_OUTPUT_ROOT`. The original study raw-count cache, sample manifest, independent identity tables, and Python dependencies used by `revision_data.py` must be present. Both cohorts are processed sequentially, one numerical thread. Completed samples are reused; partial transport files allow continuation after interruption.

Each `{dataset}/{sample}/` contains:

- `corrected_counts.npz`: SciPy CSR, **cells by genes**, with the original QC cell and unique gene-symbol order.
- `corrected_counts.mtx.gz`: Matrix Market real values, **genes by cells**.
- `genes.tsv.gz`: one `symbol` column.
- `cells.tsv.gz`: original QC metadata plus the fixed `whole_brain_broad` labels supplied to DecontX.
- `contamination.tsv.gz`: barcode, broad label, model-estimated contamination fraction, original library total, corrected estimated total.
- `audit.json`, `r_audit.json`, `input_audit.json`, `run.log`, `model_estimates.rds`, and software records.

`sample_contamination_summary.tsv` and `cell_type_contamination_summary.tsv` are descriptive summaries; cells are not treated as biological replicates. The model-estimated contamination fraction is not a physical measurement. Library-level removed-total fractions may differ from unweighted cell-level median estimates.

The corrected matrices contain fractional estimates. Do not relabel them observed UMI counts, feed them into DESeq2, or silently change detection thresholds. Normalization must use corrected totals; zero corrected totals require explicit handling. The analysis has no measured empty-droplet background, uses uncertain reference-derived broad labels, and provides a sensitivity analysis rather than proof that all ambient RNA has been removed.

The uniform maxIter=2000 numerical follow-up is specified in `../repro_decontx_convergence_amendment.md` and implemented by `../repro_decontx_extend.py`; it writes a separate directory and retains every default result. `../repro_decontx_status.py` creates a small-table convergence inventory without loading expression matrices.

Execution incident: GSM7841720 completed the official model fit (68 iterations; reported model duration 153.0552 seconds) and saved all count matrices, contamination estimates, model parameters and software records, but its audit-only script tail encountered a parsing error after the shared R source was parameterized while the interpreter was still reading it. The original script text was recovered exactly and its SHA256 matched the hash recorded before model execution. The audit was restored from the preserved native outputs without refitting the model. The original error log, recovery scripts and output checksums are retained in that sample directory. Missing end-to-end R elapsed time was left null. Subsequent invocations execute immutable per-sample R script snapshots; earlier snapshots explicitly record that they were reconstructed after execution and match the pre-run hash.

Method reference: Yang et al. *Genome Biology* 2020;21:57. https://doi.org/10.1186/s13059-020-1950-6 . Versioned method source was also retrieved from the official Bioconductor 3.23 source archive and saved in `installers/` for inspection.
