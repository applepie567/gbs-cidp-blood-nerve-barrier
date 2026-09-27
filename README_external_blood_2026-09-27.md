# Manuscript update with external blood analyses

Immune recruitment across blood cohorts and macrophage composition in inflammatory neuropathies

Updated 27 September 2026.

## Current manuscript files

- `Manuscript_refocused.docx`: revised main manuscript, five main figures and two main tables. Changed text and renumbered references are highlighted yellow.
- `Supplementary_Information_refocused.docx`: full methods and results, five supplementary tables and eleven supplementary figures, with changes highlighted yellow.
- `Figures/`: all current figures as PNG and PDF, with available SVG exports.
- `CHANGES_2026-09-27_zh.md`: explanation of the new analyses and their interpretation in Chinese.

## New evidence

GSE304872 supplies an independent cohort of four early untreated AIDP patients and three controls. Raw RNA counts were recovered for 107,656 cells, including 22,592 author-annotated monocytes. The fixed CXCL8 recruitment and complement modules had positive monocyte effects (Hedges g 0.709 and 0.882), both q 0.343. These are external estimates with directional agreement and limited precision, not statistically significant replication.

The same donors' mixed PBMC CXCL8 estimate was negative. Sorted GSE304871 cells overlap the same cohort and have preparation-dependent estimates. In the additional PRJNA1174992 PBMC reanalysis, acute CXCL8 was negative and complement positive; both q values were 0.439. Public metadata could not exclude overlap with the original GSE211225 cohort, so PRJNA1174992 was not included as an independent pooled study. Its eight paired follow-up samples showed an increased CXCL8 score, with unresolved timing and technical effects.

The targeted synthesis adds only GSE304872 monocytes to the original three cohorts and only for the two selected modules. Both resulting intervals include zero. The original seven-module synthesis is retained unchanged. The earlier GSE213455 nerve analysis is also retained; independent macrophage-specific CIDP versus CIAP replication remains unavailable.

## Data and reproducibility

Supplementary Data 1 and 2 preserve the original numerical tables and frozen v2.4.0 analysis. Supplementary Data 3 and 4 preserve the later external nerve analysis. Supplementary Data 5 supplies the new blood inputs, donor counts, sample and gene mappings, complete test tables, dated analysis plan and code. The primary scores and exact tests were independently checked in base R; the packaged analysis regenerated all five authoritative result tables exactly.

`Source_data/` provides the historical source workbook, previous nerve panel tables, new blood plot tables and an updated figure/table index. The current Figure 4 was archived Figure 3; current Figure 5 was archived Figure 4; archived Figure 5 remains Supplementary Figure S8. The historical workbook and archives retain their original labels. Follow the outer current manuscript and source index.

The original large GSE304872 Seurat object is not redistributed here. Its public download URL, size, SHA-256 and extraction scripts are supplied in Supplementary Data 5. All new statistical tables can be regenerated from the included donor pseudobulk counts without that download.

These additions have not been published as a new GitHub or Zenodo release. The DOI https://doi.org/10.5281/zenodo.22885750 refers to the earlier v2.4.0 archive and does not archive this new update. Historical documents inside the frozen archive document that earlier release; use the two outer Word files for the current text.

Public source data retain their original attribution and terms. Analysis code retains the MIT licence stated in the corresponding packages.
