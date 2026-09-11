# GSE332910 exploratory third cohort

This directory contains six newly deposited single-nucleus RNA libraries: three
Sham and three 24-hour post-reperfusion MCAO libraries. Their GEO library labels
are traceable. Supplementary Table S1 explicitly reports three biological
replicates and three RNA libraries for each group. It does not give the
individual-mouse or pooling map, so this supports author-reported biological
replicate n without independently establishing one mouse per library.
Corresponding ATAC libraries are not extra replicates.

The supplementary XLSX and DOCX were obtained from the public PMC Cloud Service
with their supplied MD5 digests verified, after discovering that PMC retired
the legacy OA API in August 2026. See
`../repro_cohort_metadata/GSE332910_supplement_download_audit.json` and the
extracted `../repro_cohort_metadata/ADVS-9999-e77547-s001.json` (Sheet1 rows 2-4).

Source article: https://doi.org/10.1002/advs.77547

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332910

The original study collected ipsilateral striatum after 60-minute MCAO followed
by 24-hour reperfusion. It used DNBelab single-nucleus libraries and MGI sequencing.
It contains astrocytes and endothelial cells, but region, assay and platform differ
from the primary scRNA-seq cohorts. Add it as a third exploratory cohort rather
than replacing the original external n=2 case group.

`../repro_third_cohort.py` downloads the six original GEO tar archives, records
SHA-256 checksums, extracts only the three expression files per archive, and runs
the original `../stream_annotate.py` QC and transparent lineage marker rule without
changing parameters. Run from any directory using an environment containing the
original analysis dependencies:

```text
python repro_third_cohort.py
```

Use `--download-only` or `--qc-only` for staged execution. Files are cached.
`archives/` preserves submitted files; `data/GSE332910/raw/` is a deterministic
filename/feature-column adapter. The one-column feature file contains submitted
mouse gene symbols. The adapter duplicates that symbol into both identifier and
symbol columns to satisfy the original reader; it does not invent Ensembl IDs.
Source matrices are gene-major sparse integer matrices, unlike the cell-major
matrices used by the previous full-CSR cache builder.

Final `library_qc_summary.tsv`, `qc_method.json` and `download_manifest.json`
record the library-level audit. `processed/GSE332910/*_cells.tsv.gz` provides the
QC metrics and annotation scores for every submitted cell/nucleus. All markers
and thresholds are in `processed/GSE332910/marker_panel.json` and `qc_method.json`.

All six archives were downloaded (561,949,676 bytes in total). There are 119,278
submitted nuclei, 116,308 passing QC, 6,430 assigned astrocytes and 1,112 assigned
endothelial cells. Each library has 75 or more assigned endothelial cells.
Of the QC-passing nuclei, 77,582 (66.70%) remain unassigned under the original
conservative marker rule. These assignments are a transferred, transparent
marker classification, not the source authors' annotations or an independent
reference-based identity validation. The `Unassigned` column in the original
QC summary also includes QC failures; use `n_qc - n_assigned` when reporting
unassigned cells among those passing QC.

After download and QC, run `../repro_third_coavailability.py`. It streams the
gene-major matrices, verifies non-duplicate cell-gene entries and assigned-cell
UMI conservation, and writes per-library raw pseudobulk, mean log1p(CP10k),
detection fractions and counts for all 15 marker classes. These are in
`processed/GSE332910/`. Full matrices are never assembled in memory.

The same script exports all 7,582 resource-covered directed astrocyte/endothelial
edges, all six library scores, gates at 5/10/20%, and the five prespecified
examples to `../../outputs/reproducibility_v2/tables/third_*`. An edge is gene
covered only if all constituent symbols occur in all six submitted feature
lists (28,483 common symbols); absent feature rows are not silently imputed.
The 10% gate retains 828 edges. Eligibility requires ligand and receptor
constituent detection of at least 10% in at least two libraries of either
condition, assessed separately for each side. Complex expression and detection
use the minimum of their constituent values. The custom score is the square
root of the product of ligand and receptor mean log1p(CP10k); its contrast is
the equal-library MCAO mean minus the equal-library Sham mean.

In the common 10%-eligible set, comparison with GSE174574 yields 141/223
same-direction edges (63.23%; Spearman rho 0.560), and comparison with GSE245386
yields 267/445 (60.00%; rho 0.547). These are descriptive transfer assessments
with different pairwise candidate universes. They do not establish a change in
reproducibility due specifically to platform or region. Table S1 reports three
biological replicates per group; individual mice per library remain unspecified.
No animal-level P values were computed. The third cohort
does not enter the original two-cohort 200-label null distribution.

For full-network LIANA, `../repro_third_liana.py --stage cache` streams just the
assigned context cells from classes with at least 30 nuclei, verifies each
retained cell's raw UMI and detected-gene totals against QC, and writes sparse
CSR counts plus the original raw column index, qualified barcode and gene order
to `context_cache/`. All six caches are complete: 4,623-8,988 context nuclei per
library, with CSR arrays of 42-134 MB. Observed cache-building peak working set
was below 475 MB. `--stage liana` and `--stage summary` are separate commands;
the first preserves native and explicitly diagnostic unique-column ranks on
the same complete cell-type network, and the latter compares all metrics on
identical complete-case candidate sets. These remain descriptive analyses.

All six LIANA runs and the summary are now complete. Their 224,004 full-network
rows and 4,194 eligible target rows, native and diagnostic ranks, source hashes,
and per-library checks are preserved in `../repro_third_liana/`. Every native
rank was reconstructed successfully, normalized stored values were strictly
positive, and the maximum observed inference working set was 1.045 GB.
The same 99 target candidates are complete across all 17 libraries. On that
fixed set, custom/native/diagnostic rho values are 0.551/0.219/0.460 versus
GSE174574 and 0.576/0.287/0.392 versus GSE245386. Same-direction counts are
53/48/52 and 58/48/48 out of 99, respectively. Native and diagnostic ranks were
computed on the full per-library networks before selecting these target edges.

`../repro_third_original_pair.py` reconstructs the original cohort-pair reference
from the same cached library scores on those same 99 targets: custom/native/
diagnostic yield 80/65/70 same-direction effects and rho 0.767/0.350/0.571.
This reference is exported separately and does not modify the previously frozen
summary. Figure 5 displays only the two third-cohort comparisons on the fixed
99 candidates; all six metrics and pairwise universes are retained in source
tables. `../repro_third_figure.py` requires inspection of its current color and
grayscale PNG before vector export. Its completed PDF, SVG and PNG use a 7.2 by
5.3 inch layout, Arial, and 7-9 pt type; PDF fonts are embedded. The third cohort
has only raw/primary results, without a reference-supported or decontaminated
version. Bilingual manuscript text is in `../repro_third_liana_narrative.json`.

Raw archive/matrix files are intentionally not required to be embedded in the code
release. Preserve the download URLs and checksums, reproducible adapter, gene
mapping, library manifest, counts and complete downstream results in the release.
