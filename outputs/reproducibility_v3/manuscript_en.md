# Expression-level concordance and implementation sensitivity of ligand–receptor priority rankings in mouse cerebral ischemia

Sihuan Zhu

Anhui University of Chinese Medicine

Correspondence: 15856965750@163.com

[ORCID 0009-0005-9705-0585](https://orcid.org/0009-0005-9705-0585)

## Abstract

Cross-cohort expression agreement can coexist with sensitivity of ligand–receptor priority rankings to software implementation and ranking context. We compared animal-level astrocyte–endothelial effects in mouse cerebral ischemia using matched candidates, dependence-preserving label permutations, controlled aggregation changes, ambient-RNA correction, a third-cohort extension and a controlled native CellChat comparison. Expression and native LIANA 1.10.0 priority concordance were 80.2% versus 56.1%; ranking each distinct score column once, matching a known upstream correction, raised priority-effect correlation to ρ=0.519; dynamically gated expression concordance was 79.9% versus a 50.0% joint-null median (exact P=0.005, exploratory-family Holm P=0.030), while fixed-metric and gap families did not reach adjusted significance. This case study supports auditing expression effects, aggregation inputs, software versions and directional denominators separately when transferring candidate priorities across cohorts.

Keywords: cerebral ischemia; ligand–receptor inference; reproducibility; astrocytes; endothelial cells; LIANA

## Introduction

Ligand–receptor inference converts cell-type expression profiles into candidate interactions and relative priorities. Cross-cohort evaluation consequently involves two related but distinct questions: whether an expression-based interaction score changes in the same direction after a perturbation, and whether the corresponding candidate gains or loses relative priority in another study. These contrasts operate on different scales and can involve different eligible interactions and network backgrounds. Agreement at one stage therefore warrants direct evaluation at subsequent stages. Distinguishing these quantities is particularly useful when a candidate list is intended to guide biological follow-up, because reproducible expression availability and reproducible prioritization represent different properties of the analysis. [1](https://doi.org/10.1038/s41467-022-30755-0)[2](https://doi.org/10.1038/s41556-024-01469-w)

Previous benchmarks have examined interaction resources, scoring methods, differential communication inference and spatial support, including reproducibility across samples or tissue sections. [3](https://doi.org/10.1093/nargab/lqaf084)[4](https://doi.org/10.1186/s13059-026-04063-5) Building on this work, an informative disease-specific assessment should compare representations on the same candidates, preserve the biological sampling unit, and identify which analytic choices account for changes in agreement. Such an assessment also requires a relevant null reference: ligand–receptor candidates share genes, so their concordance indicators cannot be treated as independent trials. Preserving each sample's complete score vector while reallocating treatment labels provides a way to retain this dependence. In parallel, explicit interventions on score representation and ranking procedures can distinguish an observed reproducibility gap from explanations that have not been directly tested.

LIANA developers merged a correction for repeated ranking of shared score columns on 9 September 2026 in [LIANA PR #261](https://github.com/scverse/liana/pull/261). We use this known correction as a controlled intervention to quantify its consequences in ischemia cohorts and the remaining dependence on ranking context.

## Methods

### Data and cell selection

GSE174574 contributed three sham and three acute middle cerebral artery occlusion (MCAO) animal libraries; GSE245386 contributed three sham and two newly generated wild-type MCAO libraries, excluding knockout samples and reused data. Samples were reconciled with GEO and the source papers, including the published correction to GSE245386. The animal represented by each original library was the comparison unit. [5](https://doi.org/10.1177/0271678X211026770)[6](https://doi.org/10.1186/s12974-023-02941-4)[7](https://doi.org/10.1186/s12974-025-03610-4)

Matrices were processed separately by library. Quality control retained 200–6,000 detected genes, at least 500 UMIs and at most 20% mitochondrial UMIs. A fixed marker rule assigned broad lineages and retained unassigned cells. In the original cohorts, an alternative selection intersected the rule labels with an agreeing CellTypist whole-mouse-brain label, reference score ≥0.5 and no Scrublet-predicted doublet. Reference scores were classifier outputs rather than calibrated identity probabilities. Full marker rules, models and diagnostic parameters are in the Reproducibility appendix. [8](https://doi.org/10.1126/science.abl5197)[9](https://doi.org/10.1016/j.cell.2021.04.021)[10](https://doi.org/10.1038/s41586-023-06812-z)[11](https://doi.org/10.1016/j.cels.2018.11.005)

GSE332910 supplied six ipsilateral striatal single-nucleus RNA libraries after sham surgery or acute MCAO. The source supplementary table reports three biological replicates and three libraries per group, without specifying mouse or pool composition within each library. The unchanged marker rule was applied, and this extension was analyzed descriptively outside the original permutation families. [12](https://doi.org/10.1002/advs.77547)

### Expression scores and relative priorities

Using a fixed mouseconsensus resource, we evaluated astrocyte-to-endothelial candidates and the reverse direction. Per-animal expression coavailability was the square root of the product of ligand and receptor mean log1p(CP10k); complexes used the minimum constituent value. Disease effects were equal-animal MCAO means minus sham means. The dynamic expression gate required all subunits to be detected in at least 10% of the relevant cells in at least two animals in either condition, assessed separately for ligand and receptor.

LIANA 1.10.0 was run per library on all eligible assigned lineages before the two target directions were extracted. CellPhoneDB, Connectome, log2FC, NATMI and SingleCellSignalR were method interfaces within this workflow. Coavailability, lr_means, expr_prod and lrscore were grouped as related expression summaries; Connectome and NATMI shared expr_prod. Relative priority was 1−magnitude_rank, where magnitude_rank is a bounded robust rank aggregation score. The interfaces therefore supplied related score representations rather than independent replications. [1](https://doi.org/10.1038/s41467-022-30755-0)[2](https://doi.org/10.1038/s41556-024-01469-w)[13](https://doi.org/10.1093/bioinformatics/btr709)

The main comparison fixed candidates with finite values for every metric in all original animal libraries, without a treatment-dependent additional gate. Directional concordance used both all candidates and, separately, those with nonzero effects in both cohorts as denominators; zero–zero pairs were not concordant. Spearman correlation included all candidates. Missing scores were not replaced by zero. Third-cohort comparisons used the same complete candidates across all original and third-cohort libraries, with rankings calculated before target restriction.

### Label randomization and multiplicity

Whole-animal score vectors were held intact while condition labels were exhaustively reassigned within each original study. The 20 discovery and 10 validation allocations yielded 200 joint configurations, including the observed labels, preserving shared-gene dependence. Each configuration recomputed concordance and correlation; the dynamic-gate analysis also reapplied expression eligibility. One-sided exact P values were upper-tail fractions of these configurations. Null percentiles describe a random-label reference interval rather than uncertainty intervals for observed concordance.

Holm correction was applied separately to family A (30 fixed-metric tests), family B (24 expression-minus-native-priority gap tests) and family C (six exploratory dynamic-gate tests). The families were specified in a revision-stage plan after earlier results had been inspected and were not preregistered. The statistics share inputs and are dependent; Holm controls family-wise error under arbitrary dependence and was retained for the original six-test family even where statistics were redundant. The gap null concerns treatment-label association, not algorithm equivalence.

A supplementary conditional diagnostic held one cohort at its observed labels while reallocating the other. Its resolution was limited by the separate allocation spaces, including a minimum validation P of 0.10. The finer joint null did not increase the number of biological replicates or establish disease association separately in both cohorts.

### Controlled aggregation and ranking context

A deterministic example and complete per-library networks were used to compare the installed aggregation with a diagnostic that ranked each distinct magnitude-score column once. We independently executed the aggregation function from the fixed merge commit of [LIANA PR #261](https://github.com/scverse/liana/pull/261) on the same saved full-network scores. This isolated the known aggregation correction; it was not a full newer-release pipeline rerun. Resource membership, detection threshold and upstream scores were held unchanged.

We then restricted aggregation to the full-network edge intersection shared by all original libraries before extracting the unchanged target set. A separate fixed-target percentile transform evaluated a simpler ranking representation. Network membership, saturation, zero effects and both directional denominators were reported. These additions were descriptive; permutation P values from the original pipeline were not transferred to them.

### Ambient RNA sensitivity

DecontX was fitted separately to all QC-passing cells in each original library using existing whole-brain reference labels as clusters and no measured empty-droplet background. Numerical convergence was checked, and the single unconverged default fit was extended under unchanged inputs and settings. All final selected fits reached the specified threshold. Model settings and the complete default-versus-extended comparison are in the Reproducibility appendix. [14](https://doi.org/10.1186/s13059-020-1950-6)

Corrected fractional counts were normalized by corrected full-gene totals while cell selections remained fixed. Detection greater than zero was the main corrected rule; estimated count ≥1 was a separate eligibility sensitivity. Original and corrected scores were also compared on identical complete candidates. These estimates were not rounded into observed UMIs or supplied to DESeq2. Count subtraction and retained eligibility remain conditional on model and cluster assumptions.

### Independent CellChat comparison

We added a controlled-resource analysis with the native CellChat R package (2.2.0.9001) using the same 11 raw-count libraries and frozen primary cell labels. The catalogue contained 1,548 shared, gene-covered interaction definitions with native complexes and cofactors. Rankings were calculated on 1,222 source–target–interaction entries with complete archived LIANA scores in every library. We compared native CellChat inferred strength, an analyst-derived strength percentile, expression coavailability and LIANA RRA recomputed once per distinct score column on this fixed network. These additional cohort contrasts were descriptive; the original permutation P values were not applied to them. [15](https://doi.org/10.1038/s41596-024-01045-4)[16](https://doi.org/10.1038/s41467-021-21246-9)

## Results

### Expression concordance against a dependence-preserving null baseline

The original cohorts retained 115,355 QC-passing cells. Under the original 10% group-dependent expression gate, 258 of 323 common directed candidates had concordant disease-effect signs (79.9%; Spearman ρ=0.743). Reapplying the gate under all 200 joint animal-label allocations gave a median concordance of 50.0% and a 95% null-reference interval of 24.5%–75.5%. Reference-supported singlet selection yielded 247/326 concordant effects (75.8%; ρ=0.703), against a null interval of 26.7%–73.3%. All six dynamically gated statistics had exact P=0.005 and Holm-adjusted P=0.030 within exploratory family C (Figure 2D). For each statistic, the observed allocation was the unique upper-tail maximum among all 200 allocations, so P=1/200 was the exact resolution limit. Within each selection, the two concordance fractions were identical across allocations because all effects were nonzero; the six-test correction was retained.

The supplementary conditional diagnostic retained the other cohort's observed labels. For both selections, dynamically gated concordance and correlation gave P=0.05 when only discovery labels were reassigned and P=0.10 when only validation labels were reassigned. These unadjusted diagnostic results distinguish the 200-allocation joint baseline from evidence attributable to each cohort separately.

### Matched candidates expose the native priority gap

The all-animal, all-metric intersection contained 187 candidates under the original selection and 174 under reference-supported singlet selection. Expression coavailability showed concordance of 150/187 (80.2%) and 127/174 (73.0%), with ρ=0.725 and 0.695. Native LIANA 1.10.0 priority effects gave 105/187 (56.1%) and 95/174 (54.6%), with ρ=0.228 and 0.102. Restricting the denominator to nonzero effects in both cohorts increased the native fractions to 105/177 (59.3%) and 95/166 (57.2%; Figure 2A–C).

Exact coavailability concordance P values were 0.010/0.020 and correlation P values were 0.005/0.005; the corresponding native-priority values were 0.145/0.225 and 0.125/0.270. The smallest Holm-adjusted values were 0.15 in fixed-metric family A and 0.12 in gap-statistic family B. Thus, the observed differences were descriptive contrasts, with family C evidence retained under its own denominator and multiplicity definition.

### A specific aggregation behavior accounts for part of the gap

A four-row deterministic example showed that the installed method loop ranked the shared expr_prod column twice, reversing its rank direction relative to the other magnitude columns. Complete-network reruns reproduced the saved target scores and native consensus in all 22 library–selection combinations. Ranking each distinct magnitude column once on these unchanged networks increased cross-cohort priority correlation to 0.519 and 0.514. Concordance became 122/187 (65.2%) and 115/174 (66.1%), or 122/168 (72.6%) and 115/158 (72.8%) among nonzero effects. These are explicitly labeled unique-column diagnostic outputs (Figure 3).

[LIANA PR #261](https://github.com/scverse/liana/pull/261) had already corrected this repeated-ranking behavior on 9 September 2026. Independently executing its fixed-commit aggregation function reproduced the unique-column diagnostic exactly on all 22 raw complete networks and the 187/174 targets (maximum absolute difference 0). This links the diagnostic to a known upstream correction while isolating the aggregation step from a complete newer release.

### Fixing the ranking universe increases correlation and zero effects

Under the original selection, complete networks ranged from 17,267–33,735 eligible edges in discovery libraries to 57,934–77,691 in validation libraries. The all-library global intersection retained 7,675 and 7,109 edges under the two selections, each covering 81 cell-type directions and containing every fixed target candidate. Unique-column RRA within these fixed universes yielded ρ=0.645/0.555. However, candidates with a zero effect in at least one cohort increased to 64/59: concordance was 100/187 (53.5%) and 84/174 (48.3%) overall, versus 100/123 (81.3%) and 84/115 (73.0%) among nonzero effects. These zero effects traced to RRA scores of 1 in every library of the relevant cohort.

A simpler fixed-target percentile transform retained expression-effect correlations of 0.609–0.697. Complete-network expr_prod and lrscore ranks were nearly redundant (ρ>0.9999999998), with small residual floating-point and tied-rank differences.

### Corrected myeloid-marker estimates and eligibility sensitivity

All eleven final DecontX fits met the specified numerical threshold. The pooled cell-level median estimated contamination was 6.13%; these cell-level summaries were descriptive. Complete fit diagnostics and the single extended fit are reported in the Reproducibility appendix.

In primary-rule target cells, mean estimated counts per cell for C1qa, Lyz2 and Tyrobp were lower than raw values in every library. Each of the six cell-type–gene combinations was summarized across eleven libraries; median library-level reductions ranged from 64.5% to 90.5% (Figure 4A,B). These changes follow model-based subtraction and do not independently validate an ambient-RNA origin. Under the >0 detection rule, corrected common sets contained 323/326 candidates. Primary and reference-supported concordance was 238/323 (73.7%; ρ=0.642) and 230/326 (70.6%; ρ=0.578). Applying an additional, stricter estimated-count ≥1 detection threshold with the same 10% group eligibility rule yielded 37/45 common candidates, with concordance of 31/37 (83.8%) and 36/45 (80.0%) and ρ=0.797/0.801 (Figure 4E,F). Within these two smaller matched candidate sets, expression-effect correlations changed from 0.795 to 0.797 and from 0.798 to 0.801 from raw to corrected data. Higher agreement in the smaller sets therefore does not establish improvement across the original candidates.

To compare the same targets, candidates complete across raw/corrected data, all eleven animals and all six metrics were fixed at 187 under primary selection and 174 under reference-supported selection. Under primary selection, raw-to-corrected correlations were 0.725 to 0.641 for coavailability, 0.228 to 0.205 for native priority, and 0.519 to 0.442 for the unique-column diagnostic. Under reference-supported selection, the corresponding changes were 0.695 to 0.660, 0.102 to 0.126, and 0.514 to 0.430 (Figure 4C,D). Animal-level scores, all-candidate and nonzero-effect denominators, and every metric are retained; raw-data permutation P values do not apply to these descriptive sensitivity comparisons. Timp3–Kdr and Ptn–Ptprz1 retained eligibility under the ≥1 threshold in both cohorts and selections, with expression effects retaining their raw-data directions. Spp1–Itga5/Itgb1 did not meet that eligibility rule in any of the four cohort–selection combinations. Plat–Lrp1 did not meet it in primary-rule discovery cells but retained eligibility in the other three combinations. Failing ≥1 eligibility does not establish absence of expression or communication, nor can the model alone determine the cellular origin of these transcripts.

### Transfer to a third cohort

GSE332910 reported three biological replicates and three RNA libraries per group. Of 119,278 submitted nuclei, 116,308 passed QC; the unchanged marker rule assigned 6,430 astrocyte and 1,112 endothelial nuclei, while 66.7% remained unassigned. Every library contained at least 75 assigned endothelial nuclei. Among candidates passing the 10% gate, third-cohort effects agreed with discovery for 141/223 candidates (63.2%; ρ=0.560) and with validation for 267/445 (60.0%; ρ=0.547). These descriptive comparisons used separate assessable intersections.

Among the same 99 target candidates complete across all 17 libraries, the original cohort pair showed 80/99 concordant custom effects (80.8%; rho=0.767), compared with 65/99 for native priority (65.7%; rho=0.350) and 70/99 for the unique-column diagnostic (70.7%; rho=0.571). Against the original discovery cohort, the third cohort yielded 53/99 for custom effects (53.5%; rho=0.551), 48/99 for native priority (48.5%; rho=0.219), and 52/99 for diagnostic priority (52.5%; rho=0.460). Against the original external cohort, the corresponding results were 58/99 (58.6%; rho=0.576), 48/99 (48.5%; rho=0.287), and 48/99 (48.5%; rho=0.392). Across the four expression-related scores, rho ranged from 0.514 to 0.720 and from 0.512 to 0.659 in the two third-cohort comparisons. The diagnostic increased priority-effect correlation, whereas direction agreement improved only slightly or remained unchanged. In comparisons involving the discovery cohort, three native and six diagnostic candidates had a zero disease effect in at least one cohort; neither priority metric had zero effects in the third-versus-external comparison. Direction agreement and effect correlation therefore capture distinct aspects of transfer. On a common target set, expression-direction transfer varied across cohorts, while priority relationships also depended on score implementation.

Barrier-related transcription, paired-pool sorted expression and spatial descriptions are provided solely as supplementary biological context (Supplementary analyses and Figures S3–S5).

### A controlled independent-framework extension

Exact native interaction definitions retained 32 of the original 187 directed candidates (29 distinct ligand–receptor definitions; 8 astrocyte→endothelial and 24 reverse-direction entries). On these same 32 entries, coavailability was concordant for 23/32 (ρ=0.565), native CellChat strength for 10/32 (ρ=0.549), the derived strength percentile for 16/32 (ρ=0.364), and matched-network LIANA RRA priority for 13/32 (ρ=0.451). Native strength had zero cohort effects for 19 entries in discovery and 3 externally; concordance among the 13 entries with nonzero effects in both cohorts was 10/13. The corresponding nonzero denominator for LIANA RRA was 19 (13/19 concordant). Three entries had zero native strength in all 11 libraries but nonzero percentile contrasts because the zero-strength tie block changed across libraries (Supplementary Figure S6).

Table 1. Data sources and analytical roles

| GEO | Contrast and sample unit | Sham / MCAO | QC retained | Role |
| --- | --- | --- | --- | --- |
| GSE174574 | Acute ischemia; animal libraries | 3 / 3 | 58,333 cells | Primary cohort |
| GSE245386 | Wild type; animal libraries | 3 / 2 | 57,022 cells | Original validation |
| GSE332910 | Ipsilateral striatum; author-reported biological replicates | 3 / 3 | 116,308 nuclei | Descriptive extension |

Table 2. Expression and priority agreement on matched targets

| Selection | Metric | N | All concordant | Nonzero concordant | ρ |
| --- | --- | --- | --- | --- | --- |
| Primary | Expression coavailability | 187 | 150/187 (80.2%) | 150/187 (80.2%) | 0.725 |
| Primary | Native LIANA 1.10.0 | 187 | 105/187 (56.1%) | 105/177 (59.3%) | 0.228 |
| Primary | Unique-column diagnostic | 187 | 122/187 (65.2%) | 122/168 (72.6%) | 0.519 |
| Primary | Fixed-global diagnostic | 187 | 100/187 (53.5%) | 100/123 (81.3%) | 0.645 |
| Reference | Expression coavailability | 174 | 127/174 (73.0%) | 127/174 (73.0%) | 0.695 |
| Reference | Native LIANA 1.10.0 | 174 | 95/174 (54.6%) | 95/166 (57.2%) | 0.102 |
| Reference | Unique-column diagnostic | 174 | 115/174 (66.1%) | 115/158 (72.8%) | 0.514 |
| Reference | Fixed-global diagnostic | 174 | 84/174 (48.3%) | 84/115 (73.0%) | 0.555 |

N denotes the fixed target-candidate count. All-candidate concordance uses N; nonzero concordance uses effects nonzero in both cohorts. Spearman correlation includes all N candidates. These representations reuse the same animals and targets and do not provide new formal between-method significance evidence.

## Discussion

This analysis identifies an implementation-sensitive component of cross-cohort ligand–receptor prioritization. Expression coavailability and native consensus priority produced different disease-effect concordance on the same candidates, and a controlled change to the aggregation procedure substantially narrowed that gap. The resulting contribution is a traceable computational case study: the score representation, eligible network and denominator can each be inspected and intervened upon. Previous evaluations have established the importance of resources, scoring choices and benchmark targets; this study adds a disease-matched assessment that connects a reproducibility pattern to a specific implementation behavior. [1](https://doi.org/10.1038/s41467-022-30755-0)[2](https://doi.org/10.1038/s41556-024-01469-w)[3](https://doi.org/10.1093/nargab/lqaf084)[4](https://doi.org/10.1186/s13059-026-04063-5)

The strongest explanatory evidence comes from preserving the complete network while changing only how shared magnitude columns are ranked. Connectome and NATMI both reference expr_prod. In the audited LIANA 1.10.0 implementation, repeated in-place descending ranking reverses this column before duplicate score names are removed. The deterministic example and full-network reproduction establish that this behavior reached the analyzed outputs. Ranking distinct columns once increased correlation under both cell selections, directly implicating this step in the observed native gap. The diagnostic exactly matches the aggregation function of the known upstream correction in [LIANA PR #261](https://github.com/scverse/liana/pull/261). It quantifies that correction's consequences while remaining distinct from a full newer-release rerun. Its improved agreement measures a computational consequence on these data, rather than establishing a biologically correct consensus.

The ranking universe supplies a second, separately controlled source of sensitivity. Fixed-global aggregation increased correlation while reducing the number of target candidates with nonzero disease effects. This combination follows from the bounded RRA score: saturation at 1 maps to zero transformed priority across an entire cohort, even when ligand and receptor expression remains measurable. RRA aggregates evidence from relative ranks and differs from a simple percentile transform. [13](https://doi.org/10.1093/bioinformatics/btr709) Reporting only conditional concordance would favor the surviving directional subset and conceal this loss of informative effects. The all-candidate fraction, nonzero fraction, correlation and zero counts together describe the intervention's consequences.

The dependence-preserving null adds a reference for expression agreement while clarifying the replication claim. Dynamic-gate concordance exceeded its joint random-label reference after correction within exploratory family C. Fixed-metric and gap families retained their separate, less conclusive adjusted results. Joint rejection concerns the configuration in which both cohorts' labels are exchangeable; it does not establish the relevant disease association separately in both cohorts. The conditional validation tail remains limited to increments of 0.10. Moving to an aggregate statistic therefore permits a finer joint null comparison while leaving the underlying small-sample design intact.

The third cohort further limits transfer. On the same 99 candidates, expression signs agreed for 80/99 between the original cohorts but for 53/99 and 58/99 in comparisons involving the third study. Tissue region, nuclear versus cellular sampling, platform and study identity change together, so their individual contributions cannot be identified here. The extension provides a descriptive boundary for this case rather than a causal explanation of cohort differences.

Ambient-count estimation adds a distinct source of sensitivity. The corrected analysis changed correlations on matched candidates, while an estimated-count threshold of at least one selected a much smaller subset. This shows why count estimation and candidate eligibility should be retained as separate outputs: an apparent gain in agreement after thresholding can reflect a change in the assessed candidates.

Native inferred strength and its derived percentile describe different features of a candidate. In this controlled extension, zero strengths reduced the number of nonzero cohort contrasts, while the percentile position of the zero-strength tie block shifted with the network background. Together with the resource attrition, these results support preserving native strengths, the exact ranking background and zero-block diagnostics when comparing candidate priorities across cohorts within an explicit resource and cell-selection scope.

Several limitations delimit interpretation. The original validation ischemia group contains two animals; the third study reports biological replication but leaves each library's mouse or pool composition unspecified. Tissue region, single-cell versus single-nucleus sampling, platform and annotation transfer accompany its more moderate expression agreement. Fixed-global reranking retains upstream expression context and therefore does not isolate cell composition. Expression-oriented metrics are also related: shared expr_prod and near-monotonic lrscore ranks provide fewer independent dimensions than the number of method interfaces suggests. Cell identity and contamination sensitivity require their own audits. Functional ligand release, receptor engagement and barrier physiology remain unmeasured.

### A practical reporting checklist

1. Preserve the biological unit: retain per-animal expression scores, condition labels and each library’s mouse or pool composition.

2. Record cell selection: publish marker panels, reference-label mappings, thresholds, doublet decisions and retained cell counts.

3. Export candidate eligibility: retain per-subunit detection, missingness and the candidate set used in each comparison and permutation.

4. Pin software and resources: record package versions, source commits, interaction resources and explicit parameters; rerun version-sensitive checks.

5. Preserve aggregation inputs: retain the complete network, distinct magnitude columns and their ranking directions before consensus calculation.

6. Report ranking context: list network membership, tie rules and saturation counts, and distinguish changing backgrounds from fixed intersections.

7. State directional denominators: report all-candidate and nonzero-effect agreement, correlation and zero-effect counts alongside all prespecified test families.

## Conclusion

The original two cohorts showed higher expression-effect concordance than native LIANA 1.10.0 priority agreement. Controlled aggregation checks reproduced a known upstream correction and explained part of this gap; ranking context and ambient-count estimates contributed further sensitivity. On the same 99 candidates, expression agreement was 80/99 between the original cohorts and 53/99 or 58/99 with the third cohort. This case supports auditing expression effects, priority effects and candidate eligibility separately, with explicit software versions and directional denominators. Broader transfer and functional interpretation require independent evidence.

## Data availability

Public source data are available from GEO: [GSE174574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174574), [GSE245386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245386), [GSE332910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332910), [GSE163752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163752), [GSE233814](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233814).

## Code availability

[https://github.com/Piggy0012/ischemia-lr-reproducibility](https://github.com/Piggy0012/ischemia-lr-reproducibility); [version DOI 10.5281/zenodo.22705840](https://doi.org/10.5281/zenodo.22705840); analysis commit `00067af9ab495f79385403141e06471b079afed6`.

The release comprises analysis scripts, versioned environment files, cell-summary counts, complete result tables and figure source data. The Reproducibility appendix and repository run instructions distinguish the full LIANA 1.10.0 workflow from replay of the aggregation function at upstream commit `d4211373692e7b9c10210488ccb1efe06452b097` on saved complete networks; the latter does not represent a full fixed-commit pipeline rerun. The file index maps each named input and output to its repository path.

## Funding

This study was supported by the Research Funds of Center for Xin'an Medicine and Modernization of Traditional Chinese Medicine of IHM (Grant/Award Number: 2023CXMMTCM002); the Anhui Provincial Clinical Translational Special Project, Research on Integrated Evaluation Model and Key Technology of "Disease-Evidence-Procedure-Effect" of Hepatolenticular Degeneration by Combination of Traditional Chinese and Western Medicine and Prevention (Grant/Award Number: 202204295107020047); and the Key Project of Clinical Research Fund of Anhui University of Chinese Medicine (Grant/Award Number: 2021sfylc08).

## Competing interests

The author declares no competing interests.

## Author contributions

Sihuan Zhu: conceptualization, research direction, critical revision, approval of the final manuscript and responsibility for the work. AI-assisted computational and writing activities are disclosed below.

## Declaration of AI use

OpenAI Codex was used to assist study planning, public-data retrieval, code development and execution, statistical and figure checks, and drafting, translation and revision of the manuscript. Numerical results were calculated by the archived analysis code from the stated source data. Sihuan Zhu reviewed and approved the manuscript and takes responsibility for its content.

## Acknowledgments

The author acknowledges the investigators who released the source datasets and the developers of the analytical software, including the LIANA maintainers for the aggregation correction documented in [LIANA PR #261](https://github.com/scverse/liana/pull/261).

## Ethics statement

This study reanalyzed public animal data and performed no new animal or human experiments. Source approvals are reported in the original studies.

## Figure legends

Figure 1. Animal libraries and uncorrected marker expression. A,B: Primary-rule astrocyte and endothelial counts; each point is one animal library and the short black line is the group mean. Discovery has 3 Sham/3 MCAO animals and validation has 3/2. C,D: Raw-count detection (>0) is calculated within each animal and averaged with equal animal weights. White separators distinguish target-lineage and myeloid markers. These are pre-correction identity audits; ambient sensitivity is shown in Figure 4.

Figure 2. Expression and native-priority effects on matched candidates and animal-label null references. A,B: MCAO-minus-Sham mean effects for the same 187 all-animal-complete primary-rule directed candidates; blue indicates nonzero concordant effects and orange denotes other outcomes. C: Observed all-candidate concordance for five metrics under primary and reference-supported selection, using 187/174 candidates. Reference segments denote the 2.5th–97.5th percentiles across 200 joint label allocations, with null medians also shown; these are not effect confidence intervals. D: The original 10% group gate is recomputed in each allocation. Observed 258/323 concordance (79.9%) is compared with a 50.0% null median; exact P=0.005 and Holm-adjusted P=0.030 refer to exploratory family C (six tests). Adjusted fixed-set results are reported separately. Whole-animal score vectors preserve shared ligand/receptor dependence. The conditional validation diagnostic still has only ten label allocations.

Figure 3. Aggregation implementation and ranking-universe sensitivity. Columns represent primary and reference-supported singlet selection, retaining the same 187/174 target candidates throughout. Native v1.10.0 is the installed LIANA consensus. Unique-column RRA ranks each distinct magnitude-score column once on the unchanged full network. The fixed-global variant additionally restricts ranking to the 7,675/7,109 complete full-network edges shared by all 11 libraries. A,B: Spearman correlation of disease mean effects. C,D: Filled points use all fixed candidates; open points use candidates with nonzero effects in both cohorts. E,F: Counts of nonzero and at-least-one-zero disease effects, which do not indicate missing expression. Absolute effects ≤10⁻¹² are treated as numerical zero. This descriptive intervention adds no P values; the bounded RRA score differs from a simple percentile rank. The unique-column diagnostic exactly matched the isolated aggregation function of the known upstream correction in [LIANA PR #261](https://github.com/scverse/liana/pull/261) on the 22 raw networks; this is not a complete newer-release rerun.

Figure 4. Ambient-RNA sensitivity. A,B: Log1p-transformed library-level mean counts per cell for C1qa, Lyz2 and Tyrobp in primary-rule astrocytes and endothelium. Thin lines connect raw and DecontX-estimated values from the same animal library (11 libraries). C,D: Primary and reference-supported selections compare cross-cohort effect correlations on the same candidates complete for raw/corrected data, every animal and all six metrics. Blue circles denote raw values and orange triangles corrected estimates. Native and Unique denote native and unique-column diagnostic priorities. E,F: Common candidate counts passing the 10% group gate under raw >0, corrected >0 and corrected ≥1 detection definitions; these eligible sets may differ. The ≥1 rule is an additional sensitivity threshold, and failure does not establish absence of expression. No new inferential tests are shown. Estimated counts, detection thresholds and numerical-convergence source selection are documented separately.

Figure 5. Effect transfer to an acute striatal single-nucleus cohort. A,B: assigned astrocyte and endothelial nuclei in each GSE332910 library under the original primary rule. Each point denotes one library; S1–S3 and M1–M3 identify the submitted Sham and 24-h MCAO replicates. Source Table S1 reports three biological replicates and three libraries per group; mice per library were unspecified. C,D: third-cohort comparisons with GSE174574 and GSE245386 on the same 99 target candidates complete across all 17 libraries. Dots show Spearman correlation of MCAO-minus-Sham effects for custom coavailability, native RRA priority and the unique-column RRA diagnostic. Labels report concordant nonzero effects/all 99 candidates and candidates with a zero effect in either cohort (absolute effect <=1e-12). Priorities were computed on each full library network before target restriction. Original-pair comparisons on the same 99 candidates are provided in the supplementary table third_rank_diagnostic_original_pair_on99.tsv. The third cohort used raw counts and primary selection only; all estimates are descriptive.

## References

1. Daniel Dimitrov, Dénes Türei, Martin Garrido-Rodriguez, Paul L. Burmedi, James S. Nagai, Charlotte Boys, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nature Communications. 2022;13:3224. [doi:10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0)

2. Daniel Dimitrov, Philipp Sven Lars Schäfer, Elias Farr, Pablo Rodriguez-Mier, Sebastian Lobentanzer, Pau Badia-i-Mompel, et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. Nature Cell Biology. 2024;26:1613-1622. [doi:10.1038/s41556-024-01469-w](https://doi.org/10.1038/s41556-024-01469-w)

3. Giulia Cesaro, Giacomo Baruzzo, Gaia Tussardi, Barbara Di Camillo. Differential cellular communication inference framework for large-scale single-cell RNA-sequencing data. NAR Genomics and Bioinformatics. 2025;7:lqaf084. [doi:10.1093/nargab/lqaf084](https://doi.org/10.1093/nargab/lqaf084)

4. Li-Ting Ku, Vincent Bernard, Jimin Min, Ying Yuan, Eugene Jon Koay, Anirban Maitra, et al. Benchmarking tools for deciphering cellular crosstalk in spatially-resolved transcriptomics. Genome Biology. 2026;27:163. [doi:10.1186/s13059-026-04063-5](https://doi.org/10.1186/s13059-026-04063-5)

5. Kai Zheng, Lingmin Lin, Wei Jiang, Lin Chen, Xiyue Zhang, Qian Zhang, et al. Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke. Journal of Cerebral Blood Flow & Metabolism. 2022;42:56-73. [doi:10.1177/0271678X211026770](https://doi.org/10.1177/0271678X211026770)

6. Zhaohui Ruan, Guosheng Cao, Yisong Qian, Longsheng Fu, Jinfang Hu, Tiantian Xu, et al. Single-cell RNA sequencing unveils Lrg1's role in cerebral ischemia‒reperfusion injury by modulating various cells. Journal of Neuroinflammation. 2023;20:285. [doi:10.1186/s12974-023-02941-4](https://doi.org/10.1186/s12974-023-02941-4)

7. Zhaohui Ruan, Guosheng Cao, Yisong Qian, Longsheng Fu, Jinfang Hu, Tiantian Xu, et al. Correction: Single-cell RNA sequencing unveils Lrg1’s role in cerebral ischemia‒reperfusion injury by modulating various cells. Journal of Neuroinflammation. 2025;22:269. [doi:10.1186/s12974-025-03610-4](https://doi.org/10.1186/s12974-025-03610-4)

8. C. Domínguez Conde, C. Xu, L. B. Jarvis, D. B. Rainbow, S. B. Wells, T. Gomes, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science. 2022;376:eabl5197. [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)

9. Zizhen Yao, Cindy T.J. van Velthoven, Thuc Nghi Nguyen, Jeff Goldy, Adriana E. Sedeno-Cortes, Fahimeh Baftizadeh, et al. A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell. 2021;184:3222-3241.e26. [doi:10.1016/j.cell.2021.04.021](https://doi.org/10.1016/j.cell.2021.04.021)

10. Zizhen Yao, Cindy T. J. van Velthoven, Michael Kunst, Meng Zhang, Delissa McMillen, Changkyu Lee, et al. A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature. 2023;624:317-332. [doi:10.1038/s41586-023-06812-z](https://doi.org/10.1038/s41586-023-06812-z)

11. Samuel L. Wolock, Romain Lopez, Allon M. Klein. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Systems. 2019;8:281-291.e9. [doi:10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)

12. Ruolin Zhang, Chang Liu, Yuneng Zhou, Kaichen Zhao, Muyang Li, Bingcheng Cai, et al. Nrsn1–Smarcc1 Coupling Regulates Neural Stem Cell Differentiation and Chronic‐phase Recovery After Ischemic Stroke. Advanced Science. 2026;e77547. [doi:10.1002/advs.77547](https://doi.org/10.1002/advs.77547)

13. Raivo Kolde, Sven Laur, Priit Adler, Jaak Vilo. Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics. 2012;28:573-580. [doi:10.1093/bioinformatics/btr709](https://doi.org/10.1093/bioinformatics/btr709)

14. Shiyi Yang, Sean E. Corbett, Yusuke Koga, Zhe Wang, W Evan Johnson, Masanao Yajima, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biology. 2020;21:57. [doi:10.1186/s13059-020-1950-6](https://doi.org/10.1186/s13059-020-1950-6)

15. Suoqin Jin, Maksim V. Plikus, Qing Nie. CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nature Protocols. 2025;20:180–219. [doi:10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)

16. Suoqin Jin, Christian F. Guerrero-Juarez, Lihua Zhang, Ivan Chang, Raul Ramos, Chen-Hsiang Kuan, et al. Inference and analysis of cell-cell communication using CellChat. Nature Communications. 2021;12:1088. [doi:10.1038/s41467-021-21246-9](https://doi.org/10.1038/s41467-021-21246-9)
