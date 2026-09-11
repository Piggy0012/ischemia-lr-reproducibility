# Running and checking the analysis

Run from the repository root with Python 3.12.14 and the [frozen environment](ENVIRONMENT.md). No command below publishes a repository, creates a DOI or submits a manuscript.

```text
python work/repro_v3_environment/verify_environment.py --release-root .
python tools/update_file_index.py --verify
python work/repro_validate.py --require-complete
python work/repro_run_cached_pipeline.py
```

The last command prints a cache-only plan. For a complete cached reconstruction in a writable working copy:

```text
python work/repro_run_cached_pipeline.py --execute --include-completed-optional --require-complete
```

All 12 stages have previously passed in both the original and relocated analysis directories, with 2,749 scientific/source checks and no scientific or protected-source differences. Their [compact evidence](outputs/reproducibility_v2/validation/relocated_cached_pipeline/) is retained. Raw/corrected matrices, model fitting and network access are blocked by the cached worker. Reruns can change gzip timestamps and provenance paths; the saved semantic comparison distinguishes these from changed scientific values. A clean conda installation or complete R lock restore was not part of those runs.

To replay the original native summary and the isolated upstream aggregation correction:

```text
python work/repro_run_cached_pipeline.py --stages raw_rank_summary upstream --execute --require-complete
```

Full cell-level inference, DecontX fitting and any independent framework are separate upstream operations requiring their own data and completion audits. Do not treat the guarded cached command as a model fit. Detailed source-scope and full-inference commands are in the [environment report](work/repro_v3_environment/report.md).

## Figures

[FIGURE_SOURCE_INDEX.tsv](FIGURE_SOURCE_INDEX.tsv) identifies actual plot inputs, exported source data, scripts, reused supplementary-figure provenance and all PNG/PDF/SVG exports for the eleven figures. Companion analyses that were not plotted have separate roles in the index. Figure scripts may write output files; inspect regenerated PNGs before exporting vectors. Optional SciPilot profiling paths in scripts refer to authoring-time tooling and are not required to read the archived figures or run the default cached tables.

## Native CellChat extension

Follow [the extension run guide](work/repro_v3_cellchat/README.md) and [its environment instructions](work/repro_v3_cellchat_environment/README.md). Full native inference was run separately for eleven libraries; it is not part of the original twelve-stage cache-only pipeline. All normalized sparse inputs and native RDS/full-network results are archived. The final numerical exports preserve exact RDS precision.

With the recorded Python environment, the following command recomputes descriptive same-universe summaries from the archived outputs without fitting a model:

```text
python work/repro_v3_cellchat_summary.py
```

Run it in a writable copy: regenerated compressed tables and timestamps can change file hashes even when scientific values agree. It reports 32 exact target entries and 1,222 network entries. It computes no new animal-level significance test. The original A/B/C tables remain under their original paths.

## Manuscript generation

Use the [separate document dependency pins](work/repro_v3_document_requirements.lock.txt) and [recorded rendering environment](work/repro_v3_document_environment.json), then run:

```text
python work/repro_v3_prepare.py
python work/repro_v3_build_docx.py
```

These commands assemble four bilingual main/supplement DOCX files under `outputs/reproducibility_v3/`, preserving the eleven reviewed figures. They read verified public identifiers from `outputs/reproducibility_v3/archive_record.json` when available and otherwise retain conspicuous draft fields. Successful document generation is followed by visual review; the original review used LibreOffice and Poppler with the specified Windows fonts. Scientific inference and manuscript formatting use separate dependency records.

## Append a completed supplement

`tools/sync_supplement.py` accepts a JSON object with a `files` array. Each item specifies `source` (an exact local file), `destination` (its intended repository-relative path) and `sha256`. Replacement of an existing different file additionally requires `expected_existing_sha256`. It prints a dry run unless `--execute` is added. It excludes runtimes, account files, large archives and files at least 100 MiB; it creates no Git commit and contacts no remote.

```text
python tools/sync_supplement.py reviewed_supplement_manifest.json
python tools/sync_supplement.py reviewed_supplement_manifest.json --execute
python tools/update_file_index.py --verify
```

Synchronize only completed code, outputs and locked dependencies. Keep publication-action drafts and private CLI state outside this directory. If a plotted input or export changes, update the figure index and perform visual review before refreshing the overall file index. Commit boundaries and publication are handled separately by the author.

Git must preserve scientific bytes: `.gitattributes` disables text conversion globally. Set `core.autocrlf=false` when initializing the repository and verify file hashes after staging; no Git initialization has been performed by the preparation scripts.
