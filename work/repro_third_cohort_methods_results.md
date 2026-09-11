# Third-cohort text for integration

## Source and methods

We extended the analysis to GSE332910, a single-nucleus RNA-sequencing study of
the ipsilateral mouse striatum after 60-min middle cerebral artery occlusion
followed by 24-h reperfusion. Six deposited RNA libraries were selected: three
Sham and three 24-h libraries. The source study's Supplementary Table S1 reports
three biological replicates and three RNA libraries per condition; individual
mouse or pool membership was not supplied. We retained these library-level
identifiers throughout the exploratory analysis. The original study used
DNBelab single-nucleus chemistry and MGI sequencing, extending the region and
assay of the two primary scRNA-seq cohorts. [Source study](https://doi.org/10.1002/advs.77547),
[Table S1](https://pmc-oa-opendata.s3.amazonaws.com/PMC13542395.1/ADVS-9999-e77547-s001.xlsx).

All six submitted integer-count matrices were downloaded with checksums, and
the original QC and marker-based primary selection rules were applied without
retuning. Source gene symbols were preserved; ligand-receptor coverage required
every constituent symbol to occur in all six feature lists. Custom coavailability
was calculated per library as the square root of the product of ligand and
receptor mean log1p(CP10k), using the minimum constituent expression for complexes.
Effects were the equal-library MCAO mean minus the Sham mean. The original
10% detection gate was applied separately to ligand and receptor constituents,
requiring at least two libraries of either condition to pass. Cross-cohort
comparisons used common eligible edges and were reported descriptively without
additional P values.

## Completed QC and expression results

Of 119,278 submitted nuclei, 116,308 passed QC. The original primary marker rule
assigned 6,430 astrocytes and 1,112 endothelial cells; every library contained at
least 75 assigned endothelial cells. However, 77,582 QC-passing nuclei (66.70%)
remained unassigned, indicating limited coverage of this fixed marker scheme
when transferred to the striatal nuclear dataset. These labels therefore provide
a consistent analysis rule, rather than independent confirmation of identity.
The intersection of submitted feature lists comprised 28,483 symbols and covered
7,582 directed astrocyte-endothelial resource edges. The 10% gate retained 828
edges. Among common eligible edges, coavailability effects agreed in direction
for 141/223 edges with GSE174574 (63.23%; Spearman rho 0.560) and 267/445 with
GSE245386 (60.00%; rho 0.547). The pairwise candidate sets differ, so their
percentages do not isolate a platform or regional contribution to transfer.

Timp3-Kdr, Ptn-Ptprz1 and Plat-Lrp1 met both the 10% and 20% detection gates in
the third cohort, with effects of +0.04331, -0.00406 and +0.27817, respectively.
The near-zero Ptn-Ptprz1 contrast warrants particular caution in sign-based
interpretation. Spp1-Itga5_Itgb1 and Col4a1-Itga3_Itgb1 had positive numerical
contrasts but failed the 10% gate, passing only the 5% sensitivity gate. Spp1
was detected in 5.41%, 11.70% and 2.40% of assigned astrocyte nuclei in the three
MCAO libraries, leaving only one library above 10%.

## Library-level QC

| Library | Condition | Submitted nuclei | QC-passing nuclei | Assigned astrocytes | Assigned endothelial cells |
| --- | --- | ---: | ---: | ---: | ---: |
| GSM9755192 | MCAO 24 h | 25,095 | 24,961 | 1,295 | 249 |
| GSM9755194 | MCAO 24 h | 21,218 | 21,154 | 1,188 | 248 |
| GSM9755196 | MCAO 24 h | 18,825 | 18,139 | 1,002 | 122 |
| GSM9755210 | Sham | 18,445 | 17,225 | 1,009 | 215 |
| GSM9755212 | Sham | 23,066 | 22,349 | 1,259 | 203 |
| GSM9755214 | Sham | 12,629 | 12,480 | 677 | 75 |

Counts refer to the unchanged primary marker rule, not author-assigned cell types.
