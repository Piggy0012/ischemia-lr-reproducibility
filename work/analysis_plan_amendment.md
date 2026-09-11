# Analysis amendment — 2026-09-09

Written after the first discovery/validation expression analysis and marker audit, before inspecting the new sensitivity results. This is an explicit post hoc sensitivity, not a preregistered analysis.

The marker audit detected myeloid transcripts in marker-identified astrocytes and endothelial cells in both datasets. We therefore add a conservative `low_myeloid` sensitivity: remove a target cell if (i) the sum of Ptprc, Tyrobp, Lyz2, C1qa, C1qb, C1qc, Csf1r and Aif1 counts exceeds 2% of its total UMI count and at least three of these genes are detected, or (ii) Ptprc has at least three UMIs. Keep the original analysis as primary. This filter can remove activated target cells; it does not estimate or remove ambient RNA and cannot establish the cell of origin of Spp1 or immune-associated transcripts.

The LIANA resource contains potential indirect or questionable biological relationships. Keep the complete resource-based results but use illustrative pairs only with explicit evidence review. A high coavailability score does not establish direct ligand binding, activation, cell contact or barrier permeability.

Spatial data use genuine sample-specific image coordinates joined to a public Visium barcode-to-array-coordinate lookup. Without the histological image and original region annotations, no ischemic-core/penumbra labels will be reconstructed. Spatial correlations are descriptive; no spot-level p-value will be interpreted as animal replication.

Small sample sizes limit exact permutation p-values. Inference about individual communication changes will explicitly prioritize this limitation over small parametric Welch p-values. Report concordance as a descriptive proportion, not a binomial test that treats dependent pairs as independent observations.
