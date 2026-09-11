# Cross-cohort ligand–receptor reproducibility assessment

Current article: **Expression-level concordance and implementation sensitivity of ligand–receptor priority rankings in mouse cerebral ischemia**.

Author: Sihuan Zhu (朱四欢), Anhui University of Chinese Medicine; ORCID 0009-0005-9705-0585. Correspondence is provided on the manuscript author page.

Public deposition is pending authenticated Zenodo access. Local archive metadata is not a DOI or proof of deposition.

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
