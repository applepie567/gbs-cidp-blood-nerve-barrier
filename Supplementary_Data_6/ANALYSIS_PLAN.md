# Nerve macrophage composition feasibility analysis

Analysis specified on 27 September 2026 before deriving the new whole-transcriptome reference or inspecting external marker-score associations.

## Question

Can a nerve-derived macrophage expression index track the macrophage contribution in held-out atlas donors well enough to justify an exploratory composition sensitivity analysis of the existing GSE213455 Fc result?

This is a post hoc extension of the manuscript. It does not create a new independent CIDP versus CIAP macrophage cohort. GSE213455 remains four CIDP and nine vasculitic neuropathy biopsies. Its existing 11-gene Fc outcome and standardisation remain fixed.

## Inputs and biological units

Use all 37 GSE285983 donors and all nuclei retained in the authors' metadata. Aggregate the deposited CellBender count matrices by donor and broad lineage. Retain donor diagnosis, centre, age and sex. The donor is the statistical unit. No nucleus is treated as an independent replicate.

Preserve the lineage distinctions in the original annotations but pool closely related clusters: macrophages, Schwann cells, perineurial cells, endoneurial stroma, epineurial stroma, blood endothelial cells, lymphatic endothelial cells, mural cells, T and NK cells, B cells, mast cells, granulocytes and adipocytes. Individual donor-lineage reference profiles require at least 20 nuclei. Whole-donor mixtures retain all annotated nuclei.

## Fixed index construction

Use genes with an unambiguous single-symbol mapping on GPL13369. Exclude all 11 Fc outcome genes, mitochondrial and ribosomal genes, XIST and common Y-linked markers. These exclusions reduce direct score overlap but do not make the marker index independent of macrophage activation.

Within each training set, average donor-level CPM profiles with equal donor weights separately for each lineage. Select up to 20 macrophage markers with macrophage mean CPM at least 10, expression above 1 CPM in at least half the eligible macrophage donors, and at least fourfold enrichment over the highest other-lineage mean after adding 1 CPM. Rank eligible genes by that enrichment. Require at least 10 selected genes.

For each whole-donor or whole-biopsy expression profile, rank all common unambiguous genes other than the 11 Fc genes within that sample. The macrophage index is the mean percentile rank of the selected marker genes. It is an expression index, not an absolute cell fraction.

## Feasibility gate

Leave out one entire atlas donor at a time. Rebuild marker selection without that donor and predict its index. Compare this index with the observed fraction of macrophage nuclei and the macrophage share of captured counts. Require Spearman correlation at least 0.70 with captured-count share, at least 0.60 with nucleus fraction, and positive correlations within each centre that has at least five donors. Also repeat the construction leaving out each entire centre. Require positive centre-held-out correlations with both targets in every evaluable centre.

If these conditions fail, do not use the marker index as a composition adjustment covariate. Report the feasibility result and stop the external adjustment. Do not tune marker thresholds or genes against the external Fc effect.

## External analysis if the gate passes

Freeze markers from all atlas donors. Calculate the index in the 13 eligible GSE213455 biopsies from the deposited log-scale expression matrix, taking the median of probes with exact single-symbol mappings. Fit Fc score against diagnosis and this one index. Use an HC3 sandwich interval with residual t degrees of freedom. Retain the original unadjusted result as the primary disease comparison. The new conditional model is exploratory and does not estimate expression within measured macrophages.

Report probe-mean aggregation and leave-one-donor-out coefficient ranges as sensitivity analyses. Report score correlation, leverage and group overlap. Do not select the strongest result. Do not interpret conditional attenuation as a causal mediation estimate. Cross-platform accuracy and cell-state effects remain unverified even if atlas calibration passes.

## Outputs

Supply donor-lineage aggregate counts, annotation and sample maps, the full probe mapping, marker lists for all held-out analyses, feasibility metrics, any permitted external models, a figure with source tables, and executable scripts. Add concise manuscript wording and full supplementary methods/results. Preserve the existing blood results and the frozen public archive.
