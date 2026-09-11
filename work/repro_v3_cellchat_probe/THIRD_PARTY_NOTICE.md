# CellChat source and native mouse resource

The files in this directory attributed to CellChat, including copied R source, DESCRIPTION, NAMESPACE, vignette and native CellChatDB.mouse data and their tabular exports, retain the upstream CellChat licensing. They are not relicensed under the manuscript project's MIT license.

The production comparison uses the complete native CellChat package version 2.2.0.9001 built from the official jinworks/CellChat commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`. Upstream DESCRIPTION specifies GPL-3. The full license text is preserved as `LICENSE.CellChat`, with its download provenance in `license_provenance.json`. Source URLs and SHA256 values are recorded in `source_manifest.json`; package binary provenance is recorded separately in `../repro_v3_cellchat/download_manifest.json` and `installed_package_verified.json`.

The full official GitHub source archive for this exact commit is separately preserved as a release asset, with URL, SHA256 and checks against the pinned DESCRIPTION, core R source and native mouse database in `corresponding_source_asset.json`. The archive contains the full upstream commit tree; it supplements the smaller reviewed source subset stored directly in this analysis repository. `../repro_v3_cellchat_fetch_source.py` restores the exact source asset.

The `stable_v2.1.2` subdirectory contains an earlier official version inspected during feasibility testing; it was not the package used for the production models. Its original DESCRIPTION and source manifest are retained. The source-only feasibility toy is not a substitute for the complete native package used in the eleven final models.

The repository's original analysis scripts, resulting analysis tables and provenance records should be distinguished from these third-party source and database materials. See the upstream repository for attribution and applicable redistribution terms: <https://github.com/jinworks/CellChat>.
