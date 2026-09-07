# Version 41 exploratory strengthening plan

Written 2026-09-07 before running the new sensitivity analyses. This is a transparent analysis specification, not a preregistration. The macrophage Fc receptor association from version 39 was already known and motivated these checks.

## Donor level nerve analysis

Use the archived GSE285983 donor pseudobulk tables with at least 20 nuclei per donor and cell compartment. Reproduce the original CIDP versus CIAP macrophage Fc receptor module result before extending it. Preserve the original gene standardisation across eligible donors in all diagnostic categories for the main sensitivity analyses.

Report all 12 component genes with a Benjamini–Hochberg family of 12. Also report four exploratory functional summaries: activating receptors FCGR1A, FCGR2A and FCGR3A; inhibitory receptor FCGR2B; neonatal Fc receptor FCGRT; and signalling components FCER1G, TYROBP, SYK, LYN and HCK. FCGR2C and FCGR3B remain in the original 12 gene score and in the gene table, but are not assigned to these four summaries. This grouping is an interpretive sensitivity analysis selected after the original result.

Evaluate the original module with two sided Mann–Whitney tests in all donors and male donors, donor case-comparator-only standardisation, ordinary least squares with age and center, and with age, center and sex. Use HC3 covariance with residual degrees of freedom t intervals. Adjustment cannot estimate a female CIDP effect because all CIDP donors are male. Report a two sided exact center-stratified label permutation using a center-adjusted disease coefficient and fixed case counts within each center. Its exchangeability assumption does not establish causality.

Report leave-one-donor, leave-one-center and leave-one-gene sensitivity estimates without calling these independent replication. Expand multiplicity correction across every eligible cell compartment and module for the CIDP versus CIAP contrast as a sensitivity check. Do not substitute the broader family retrospectively for the original declared family without disclosure.

## Subtype and external data feasibility

Download the original immune annotation and inspect macrophage subtype composition at the donor level. Test all reported macrophage subtypes as one family. Where raw counts can be obtained, test whether the Fc score persists after accounting for macrophage composition. Spatial measurements from the same atlas are within-study triangulation, not external replication. Check measured gene coverage, donor groups and access before defining spatial tests.

Check PRIDE PXD056286 for donor-resolved quantitative protein tables. Reanalyse only if sample identities, comparator labels and abundance values are verifiable. Raw instrument files alone are not equivalent to a ready quantitative matrix. Inspect public CIDP blood deposits and preserve small discovery sample sizes. Do not count unrelated healthy controls pooled from other studies as independent balanced validation.

Report null results and access limitations. Do not fabricate unavailable clinical variables, subject counts or experimental evidence. Public outputs contain aggregate results, code and original-source links. Donor intermediates are retained for private author use under the existing publication restriction.

## Coverage clarification before the first successful statistical run

The archived expression table contains 11 of the 12 specified Fc genes. FCGR2C is absent. The original score averages available genes, so the reproduced score and component gene family contain 11 measured genes. FCGR2C will be explicitly recorded as unavailable and will not be assigned a fabricated expression value or a statistical test. This clarification was made after the first run stopped at an input coverage error, before generating inferential results.

## Subtype specification before subtype results

The immune annotation is a refined reclustering of the original atlas. Barcode joins show that some nuclei originally called macrophages were subsequently assigned to dendritic, granulocyte or other populations, and some lack a refined annotation. Retain these categories rather than silently discarding them. Report every observed refined category within the original macrophage denominator, with one BH family. For composition adjustment use the first two principal components of donor centered log ratios, adding 0.5 to each category count, before fitting the existing score with age and center. Treat this as a descriptive sensitivity analysis that cannot distinguish confounding from disease-mediated composition changes.

Reaggregate all 20 case and comparator donor matrices. Reproduce the original broad macrophage log CPM, then analyse the original Macro1 and Macro2 groups and the refined union Macro1 through Macro18. For individual refined macrophage subtypes retain at least 20 nuclei per donor and at least three donors per diagnosis. Standardise genes across eligible CIDP and CIAP donors within each subtype. Report all eligible subtype tests as one BH family. These post hoc tests are not independent replication.

Spatial H5 coverage inspection shows that none of the 12 specified Fc genes are in the 99 gene panel. Limit this dataset to an explicit coverage and donor-count audit. Do not infer unmeasured Fc expression from marker proximity.
