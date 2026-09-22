# Focused Fc receptor follow up

Date 2026-09-21

This exploratory follow up was planned after the previous manuscript results were known. It is not a preregistered confirmatory study. Decisions below were recorded before calculating the new expression contrasts. Source file inventories and nucleus counts had already been inspected.

Provenance clarification added after source-code review. The original pseudobulk pipeline already calculated CIDP versus CTRL comparisons. This follow up reconstructs and explicitly reports that Fc comparison, adds CIAP versus CTRL, and interprets both alongside the reported CIDP versus CIAP contrast. The two reference comparisons are a focused reporting family, not two wholly unexamined tests. The 91-test family below is a sensitivity analysis of the manuscript's 89-test family plus these two reference comparisons, not a correction across every analysis ever run in the original pipeline.

## Protein evidence

Inspect the public PXD056286 Proteome Discoverer database for explicit specimen labels, treatment or diagnosis mappings, protein identities, unique peptides, confidence filters and quantitative values. Instrument files are not assumed to be independent patients. Proceed to disease comparisons only if the database supplies an unambiguous biological specimen and group mapping. Do not infer diagnoses from filenames or use the manuscript sample total to invent mappings. Verify quantitative field meaning before calculating effects. Keep all measured genes from the existing Fc receptor, endothelial and Schwann cell panels in the coverage audit. Missing proteins are unmeasured, not evidence of no change.

## Nerve reference comparisons

Retain the original 11 measured Fc genes, macrophage annotation, 20 nucleus threshold, log2(CPM + 0.5) expression, and 37 donor standardisation. Recompute the existing scores from archived donor expression and require agreement with the archived scores and CIDP versus CIAP result. Add exactly two comparisons, CIDP versus CTRL and CIAP versus CTRL. Use two sided exact Mann Whitney tests, BH across these two new comparisons, and a sensitivity family containing the previous 89 nerve tests plus these two tests. Describe the three diagnosis groups together with the original contrast as context. Estimate unadjusted mean differences with patient bootstrap percentile intervals, 20000 replicates, seed 20260921. The two new contrasts are descriptive given four graft reference donors, centre imbalance and sex imbalance. Do not interpret a nonsignificant comparison as equivalence or a normal expression level.

## Composition standardisation

Retain the original 20 nucleus eligibility threshold for each donor and refined macrophage subtype. Audit all 18 subtypes in the 9 CIDP and 11 CIAP donors. Direct standardisation across all 20 donors requires at least one eligible common subtype and nonzero common support for all included strata. If no subtype meets the existing threshold in every donor, report that direct full-cohort standardisation is not estimable at this resolution. Do not fill missing subtype expression, lower the threshold after inspection, or claim a composition adjusted intrinsic effect from incomplete support. Retain and appropriately qualify the existing composition regression and subtype tests.

## Reporting

Report results regardless of direction or significance. Add supported results and the measurement and eligibility findings to the manuscript and analysis attachments. Keep the existing public repository version and journal formatting for the later requested stages. No external upload or author contact is part of this analysis.
