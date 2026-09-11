# DecontX ambient-RNA sensitivity amendment

Written 2026-09-11 before any DecontX outputs were examined. This is an amendment to an exploratory reanalysis, not preregistration.

All 11 original QC-passing libraries in revision_sample_manifest.tsv are analysed independently, retaining all QC-passing cell types. Input is the original nonnegative integer UMI count matrix from revision_cache, with duplicate gene symbols summed by revision_data.load_sample. No corrected matrix is called raw counts.

Use the official Bioconductor decontX package in a workspace-local R installation. Supply z = the previously computed whole_brain_broad reference labels (all cells, without probability filtering), background = NULL, batch = NULL because each invocation contains one library, seed = 20260911, delta = c(10,10), estimateDelta = TRUE, maxIter = 500, convergence = 0.001, iterLogLik = 10. Other parameters remain the installed-version defaults and are saved. Default varGenes/dbscanEps do not control clustering because z is supplied. No disease-condition information is passed to the model. Keep rare labels as supplied; do not select labels or thresholds based on candidate results.

Run one sample per process, one numerical thread. Save counts, gene/barcode order, per-cell estimated contamination, input hashes, software versions, parameters, fitted estimates, sessionInfo, logs, and checks for finite nonnegative corrected values. Save corrected counts as sparse Matrix Market real-valued data plus portable CSR NPZ. They are estimates and can be fractional; do not round them to create fictitious observed counts or feed them into DESeq2. Subsequent expression-score sensitivity analyses may use these fractional estimates with their correct provenance.

The model infers contamination from filtered cells without empty droplets. These inferred quantities cannot establish measured ambient composition or complete removal; they depend on cluster identity assumptions and can remove shared native expression. Reference labels were learned in healthy mouse reference tissue and are uncertain in ischemia. Retain original analyses and report the DecontX analysis as sensitivity analysis, including candidate losses and disagreements. This task does not perform SoupX or claim agreement between two ambient-RNA methods.

Official sources consulted: https://bioconductor.org/packages/release/bioc/html/decontX.html ; https://bioconductor.org/packages/release/bioc/manuals/decontX/man/decontX.pdf ; https://cran.r-project.org/bin/windows/base/ .
