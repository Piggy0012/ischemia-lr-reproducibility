# Analysis plan version 1

Frozen before examination of expression differences, 2026-09-09. This is a timestamped local analysis plan, not a registered protocol.

## Objective
Evaluate reproducibility of astrocyte–endothelial candidate ligand–receptor expression changes following experimental cerebral ischemia and their relationship to barrier-associated endothelial transcription. Pure public-data analysis; no inference of causal signaling or functional barrier permeability from RNA alone.

## Selection and inference
- Discovery candidate GSE174574, 24-hour MCAO versus sham, three deposited libraries in each group; biological replication to be confirmed against primary methods. Alternative or additional sc/snRNA datasets only enter confirmatory inference if independent biological units and relevant target cells are traceable.
- GSE163752 is supportive cell-type expression validation (sorted populations, ipsilateral versus contralateral at 24 hours), not validation of single-cell subtype identity or physical communication. Pool structure and sample mapping require audit.
- Spatial data contribute within-section regional support; mixed spots are not individual cells. Sections from the same animal are not independent animals; a single section per condition precludes population-level between-condition significance.
- Distinguish MCAO model, reperfusion, time after insult, region, sex, age and control type. Do not merge disparate stages into a continuous biological trajectory. Hold out validation data from discovery selection.
- Primary species: mouse. Other species may provide supportive direction checks, never direct pooled expression analysis.

## Processing and transparent identity annotation
Use submitted raw UMI count matrices. Process one sample at a time under the available memory limits. Retain cells with 200–6000 detected genes, >=500 UMIs and <=20% mitochondrial UMIs; report attrition per sample. Repeat the key conclusions with <=10% mitochondrial UMIs and stringent cell identity as sensitivity analyses. The mitochondrial threshold is a sensitivity choice, not evidence that discarded ischemic cells lack biological relevance.
Identify astrocytes/endothelial cells from independently sourced lineage marker panels and evaluate off-target lineage expression. Keep identity markers separate from barrier gene-set outcomes where possible. Do not assign cell identities based on the desired disease association. Save all thresholds, ambiguous cells and annotation diagnostics. Clustering/embedding is descriptive and not an independent replicate.

## Expression and program analyses
Aggregate raw counts by biological sample and cell type (pseudobulk). Require >=30 cells for the primary sample×cell-type analysis, with >=20 and >=50 thresholds as sensitivity checks. If this removes a group down to fewer than three biological units, downgrade the analysis to exploratory rather than manufacturing replication.
Use established count-based differential expression (PyDESeq2) with the study-appropriate design, BH FDR and effect sizes. An expression filter (at least 10 raw counts in at least the smaller group's number of samples) is fixed independently of outcomes. Check dispersion, normalization, outliers and genes excluded by filtering.
Barrier-maintenance, endothelial-activation and angiogenesis programs will be defined from primary sources before examining condition differences. Report per-sample scores and individual genes; do not call a transcriptional score a direct barrier function measurement. Program associations are exploratory with small biological n, especially when driven by condition differences.

## Candidate ligand–receptor analysis
Freeze official LIANA mouseconsensus resource and checksum. Both directions (astrocyte→endothelial and endothelial→astrocyte) are eligible. Require all members of a complex to be observed; report detection fractions, expression and response at the sample level. The candidate universe is restricted by expression/detection independent of differential p-values. Predefined expression gate: each subunit in >=10% of relevant cells in >=2 biological samples in at least one group. Evaluate 5% and 20% gates as sensitivity checks.
Candidate scores summarize ligand and receptor abundance; they are not measured signaling probabilities. Independent replication is assessed for constituent expression changes, effect directions and candidate ranks as permitted by each validation dataset. Unreplicated, opposite-direction and already-established candidates remain in the report. A source used to define a score does not simultaneously validate that score.

## Outputs and reproducibility
Deliver source/sample manifest, inclusion/exclusion log, software versions, scripts, complete tested-gene and candidate tables, figure source data, analysis-derived manuscript and limitations. Distinguish retrospective data reuse from original animal experimentation. Do not invent author names, affiliations, ethics approval, funding, conflicts or experiments. The manuscript cannot be described as submission-ready without author-provided declarations and final author review.
