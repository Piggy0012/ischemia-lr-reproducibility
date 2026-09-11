# Independent CellChat controlled-resource extension

This directory contains a separate exploratory analysis using the **complete native CellChat R package**, not a LIANA wrapper. Read `../repro_v3_independent_framework_plan.md` and `../repro_v3_cellchat_precision_amendment.md` before interpreting results. Primary scope is original raw/primary libraries from GSE174574 and GSE245386, with a fixed shared 1,548-pair mouse resource. Every animal/library remains a separate model input. No new disease-level P values are produced.

The package is version 2.2.0.9001, development upstream commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`, built for Windows x86_64 R 4.6.1 by r-universe from the official jinworks repository. Package identity, 47 function definitions, complete native mouse DB equality, and the installed-package toy are verified in `installed_package_verified.json`. The independent `rlibrary` takes precedence only for these analysis processes. The former DecontX library supplies previously installed, read-only dependencies.

## Reproduction on the archived Windows R environment

Run from the project root with the archived `work/revision_env/Scripts/python.exe` and `work/r_runtime/bin/Rscript.exe` (or equivalent explicitly selected interpreters).

1. `python work/repro_v3_cellchat_fetch.py` verifies/downloads the exact binary files in `download_manifest.json`. If a rolling upstream CRAN URL has expired or changed, use the archived ZIP identified by its SHA256; do not silently substitute a newer version.
2. `Rscript --vanilla work/repro_v3_cellchat_install_verify.R` installs into this directory's isolated `rlibrary`, then verifies the full package, exact official core source and native mouse DB and runs the toy.
3. `python work/repro_v3_cellchat_run.py --dataset GSE245386 --sample GSM7841720` exports original normalized signaling input and runs the largest controlled library with a 2 GiB measured R process-tree RSS stop.
4. `python work/repro_v3_cellchat_run.py --all` runs or validates the remaining controlled libraries sequentially. Completed models are cached with provenance.
5. `Rscript --vanilla work/repro_v3_cellchat_export_lossless.R` exports all completed native networks through verified 17-significant-digit roundtrips.
6. `python work/repro_v3_cellchat_summary.py --prepare-only` fixes the common full-network keys using existing LIANA eligibility and the shared resource; `python work/repro_v3_cellchat_summary.py` compares all completed native models on those exact keys. In the original execution, the preparation occurred before inspecting new CellChat cross-cohort results.
7. `python work/repro_v3_cellchat_attrition.py` regenerates the edge-by-edge explanation of the 187 original candidates narrowing to 32 exactly represented native CellChat pairs.
8. `python work/repro_v3_cellchat_figure.py` rebuilds the S6 PNG and figure source tables. Local SciPilot authoring helpers are optional; without them a public pandas-based profiling fallback produces the identical plot. `--export` additionally writes PDF/SVG after the PNG matches its completed visual-review record.

The input exporter recreates the prior float32 logCP10k operations with all-gene UMI denominators, then promotes values losslessly for native R storage. Only signaling genes and frozen context cells are exported; no full-gene dense matrix is made. Each folder under `controlled/<GSE>/<GSM>` contains the gene/cell order, normalized sparse input, input audit, native model RDS, all-row score exports, original internal P values, zero scores, execution log and per-second memory trace.

Use `full_network_lossless.tsv.gz` for numerical analysis, not the initial `full_network.tsv.gz`. The native RDS is the authoritative full-precision model output. The later CellChat strength percentile is an analyst-derived relative rank, **not a native CellChat consensus rank**. Native cell-label permutation P values do not measure cross-animal disease effects and are not used to filter the controlled comparison.

Archive source files, full networks, target/sample scores, complete comparisons, resource/mapping/attrition tables, package/runtime versions, download hashes, and logs. Do not put unpacked `rlibrary` trees or original raw sequencing matrices into the code archive. Preserve the 42 compact binary ZIPs in a separate environment archive if durable Windows reinstallability is required; their rolling upstream URLs alone do not guarantee future availability.

## Completed comparison and release contents

All 11 models and full-precision probability/P-value exports passed their checks. The measured maximum R process-tree RSS was 1,417,273,344 bytes (1.32 GiB). On the same 32 target entries, expression coavailability had rho 0.565 and 23/32 concordant effects; native CellChat strength had rho 0.549 and 10/32, with 10/13 concordant among effects nonzero in both cohorts; the analyst-derived strength percentile had rho 0.364 and 16/32; matched-network LIANA unique-column RRA had rho 0.451 and 13/32 (13/19 among effects nonzero in both cohorts). These are descriptive results conditional on exact resource representation and archived LIANA eligibility.

Three targets have zero native CellChat strength in every library but nonzero derived percentile contrasts as the zero-strength tie block changes with the network background. This is a property of relative rankings, not evidence for altered biological communication. Full zero/tie audits and all effect directions, including failures and opposing directions, are retained in `tables`.

The extension archive index is `../../outputs/reproducibility_v3_cellchat/cellchat_archive_manifest.json` with a TSV companion. It includes 865-gene normalized sparse inputs, native RDS networks, all-row native strength/internal P exports, exact comparison inputs, scores, source tables, code, fixed resource definitions, licenses, environment metadata and independent audits. Runtime/package binary assets are deliberately listed separately for release storage. The effective environment and actual offline restoration audit are in `../repro_v3_cellchat_environment`. Upstream CellChat source and database licensing remains GPL-3 as recorded in `../repro_v3_cellchat_probe/THIRD_PARTY_NOTICE.md`; the project license does not relicense third-party content.
