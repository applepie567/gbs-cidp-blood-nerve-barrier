# Release v2.3.1 — 16 September 2026

This release synchronises the manuscript, aggregate source workbook, evidence tables and Figures 3–5 following source verification.

- Figure 3C groups CSF concentration measurements as one assay type. Protein accumulation is supported by reported concentration measurements, not by the identification count of 1402 proteins. The Ding et al. entry records 97 differential proteins among 298 identified proteins (30 higher and 67 lower).
- Figure 4A uses equal-area dots and mean donor log2(CPM + 0.5) in seven cell groups with matching summaries. Figure 5B uses mean donor detection percentages and mean donor log2(CP10k + 1), weighting nine CIDP donors equally.
- Figure 5A states the evidence type in each cell. Statistical tests with different correction families and descriptive expression localisation are not ranked on one scale. All blood pooled effects have P > 0.05. The macrophage Fc result has q = 0.0087 within eight modules and q = 0.087 across all 80 nerve module comparisons.
- The nerve interferon entry identifies OAS1 enrichment in macrophages and MX1 enrichment in arterial endothelial cells in Heming et al. Supplementary Data 3, and lower MX1 in nonmyelinating Schwann cells in polyneuropathy versus controls in Supplementary Data 9. This is not a dedicated CIDP–CIAP interferon-module test.
- Table 3, the cross-compartment map and published nerve table are synchronised. C3, C5 and C9 were reported as higher in the 9-CIDP/2-control proteomic subset; C6 was elevated in only two cases. The 52/55 biopsy deposition count remains a separate observation.
- Current Word files contain the revised English text and updated figure legends. Review comments are excluded from public copies. Historical v41/v42 files remain explicitly labelled historical.

The blood and CSF pooled estimates, nerve module estimates and genetic contrasts are retained. This release does not claim a new raw-matrix analysis or external validation. Additional plot inputs are group aggregates. The donor-level genetic input in manuscript Additional file 4 is not included in this public aggregate archive. Original database matrices remain at their cited sources.

Run `python analysis_update/analysis/run_v231_figures.py` for the revised figures. The public aggregate verification script checks the retained supported calculations. Code: MIT. Project-derived data and figures: CC BY 4.0. External sources retain their original terms.
