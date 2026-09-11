# Ischemia ligand–receptor reproducibility

**Expression-level concordance and implementation sensitivity of ligand–receptor priority rankings in mouse cerebral ischemia**

Sihuan Zhu, Anhui University of Chinese Medicine. [ORCID 0009-0005-9705-0585](https://orcid.org/0009-0005-9705-0585).

This repository contains the audited analysis code, complete result tables, cell-summary counts, figure inputs and environment records for this reproducibility study. Version 3.0.0 preserves the original v2 numerical analysis and adds an independently executed native CellChat comparison. The public code/results release and the final manuscript are separate versioned artifacts.

Start with [RUNNING.md](RUNNING.md), [ENVIRONMENT.md](ENVIRONMENT.md), the searchable [file index](FILE_INDEX.tsv), and the [eleven-figure source index](FIGURE_SOURCE_INDEX.tsv). The [v2 result guide](outputs/reproducibility_v2/README.md) describes the original statistical families; its manuscripts are retained historical drafts. Earlier tests and source drafts remain available under their original versioned paths. The v3 manuscript generators and audited bilingual text additions are included in `work/`.

The [reviewed revision-3 manuscripts](outputs/reproducibility_v3/README.md) include English and Chinese main/supplement DOCX and PDF files. They are explicitly drafts pending the version-specific Zenodo DOI; no preprint or journal submission has occurred. Their scientific analysis is pinned to commit `00067af9ab495f79385403141e06471b079afed6`, separately from subsequent editorial commits. The release and environment/source assets are prepared as a draft awaiting the author's Zenodo connection. See the [publication status](outputs/reproducibility_v3/publication_status.json) and the [verified upstream feedback](https://github.com/scverse/liana/pull/261#issuecomment-5632198130).

## Direct data entry points

- Fixed cell-identity rules: [discovery marker panel](work/processed/GSE174574/marker_panel.json), [validation marker panel](work/processed/GSE245386/marker_panel.json), [third-cohort marker panel](work/repro_third_cohort/processed/GSE332910/marker_panel.json).
- Original cohort comparison restricted to the added cohort's common 99-candidate set: [full comparator table](outputs/reproducibility_v2/tables/third_rank_diagnostic_original_pair_on99.tsv), [animal-derived effects](outputs/reproducibility_v2/tables/third_rank_diagnostic_original_pair_on99_effects.tsv.gz), [audit](outputs/reproducibility_v2/tables/third_rank_diagnostic_original_pair_on99_audit.json). This companion analysis is distinct from the plotted original-versus-third comparisons.
- Sorted data: [all 49,914 paired-pool gene effects and intervals](outputs/reproducibility_v2/tables/sorted_all_gene_descriptive_effects.tsv.gz), [fixed illustrative component effects](outputs/reproducibility_v2/tables/sorted_fixed_component_effects_ci.tsv), and [retained original tests](work/sorted_rna/diff_expression.tsv.gz).
- [Null distributions and fixed candidate data](outputs/reproducibility_v2/tables/), [figure source data](outputs/reproducibility_v2/figure_source_data/), [all-library raw/corrected LIANA network exports](work/repro_liana_rank_diagnostic/), and [third-cohort networks](work/repro_third_liana/).
- [Successful relocated cached reconstruction](outputs/reproducibility_v2/validation/relocated_cached_pipeline/) and [publication-package provenance](provenance/v2_archive/release_integrity.json).
- Independent native CellChat: [run instructions](work/repro_v3_cellchat/README.md), [complete controlled models](work/repro_v3_cellchat/controlled/), [all comparison tables](work/repro_v3_cellchat/tables/), [S6 and source data](outputs/reproducibility_v3_cellchat/), and [independent numerical review](work/repro_v3_cellchat_independent_review.md).

## Independent framework and candidate scope

CellChat 2.2.0.9001 was executed as a complete native R package on the same eleven original libraries, using a controlled shared resource. Exact native ligand/receptor definitions retained 32 of the original 187 directed target entries; the remaining 155 were resource-definition exclusions. Relative scores used the same 1,222-entry network. Native inferred strength, an analyst-derived strength percentile, expression coavailability and unique-column LIANA RRA are separate outcomes, with complete and nonzero-effect denominators reported. These descriptive results do not inherit the original A/B/C permutation P values. Native all-row scores, internal P values, zero blocks and lossless RDS exports are retained.

## Interpretation and attribution

LIANA 1.10.0 complete per-library inference and a controlled aggregation correction are different execution scopes. The correction was already merged upstream in [LIANA PR #261](https://github.com/scverse/liana/pull/261), commit `d4211373692e7b9c10210488ccb1efe06452b097`. Four aggregation functions were replayed on the same saved networks; a complete newer-version cell-level pipeline was not run. This case evaluates consequences of the known correction and remaining context sensitivity, not discovery of the software defect.

Statistical units, fixed candidate universes, both direction denominators, the joint and conditional nulls, and separate multiplicity families are preserved. Corrected fractional counts, implementation diagnostics and the third cohort do not inherit native raw P values. See the [analysis plans](work/repro_null_plan.md) and [implementation note](work/repro_upstream_correction_note.md).

## Included data and integrity

The original ZIP's 1,596 payload files and 269 original v3 environment files are retained alongside the completed CellChat extension and its own environment manifest. Native relative paths and all original scientific bytes are preserved. Repository-entry README, citation and Zenodo metadata are updated; their original versions are archived under [provenance/v2_archive](provenance/v2_archive/). The original `MANIFEST_SHA256.tsv` retains its original ZIP scope. [FILE_INDEX.tsv](FILE_INDEX.tsv) is the current repository inventory; it excludes itself and its external checksum to avoid circular hashing.

Raw GEO downloads and large corrected cell matrices are obtained separately using the preserved accession manifests and checksums. This repository includes their counts/aggregates, complete results and audits, not a duplicated runtime or raw-cell data mirror. No included file reaches GitHub's 100 MiB regular-Git threshold. The matrix asset exceeds 2 GiB as a single ZIP and needs Zenodo or smaller cohort assets; see the [size and environment report](work/repro_v3_environment/report.md).

Original analysis code uses the [MIT license](LICENSE); derived tables, figures and manuscripts use CC BY 4.0 as stated there. Third-party source, datasets and resource licenses retain their own terms. Use [CITATION.cff](CITATION.cff) and cite the original dataset publications and method developers. [Author declarations](AUTHOR_DECLARATIONS.md) record the author-supplied funding and competing-interest statements.
