# Immune gene and protein profiles across blood, cerebrospinal fluid and peripheral nerve in inflammatory neuropathies

Release **v2.3.1** contains the source-verified revision dated **16 September 2026**.

[GitHub release](https://github.com/applepie567/gbs-cidp-blood-nerve-barrier/releases/tag/v2.3.1) · [Zenodo version series](https://doi.org/10.5281/zenodo.22226672)

| Material | Current location |
|---|---|
| Revised manuscript with five embedded figures and three tables | [Word](GBS_CIDP_Manuscript_v2.3.1.docx) · [PDF](GBS_CIDP_Manuscript_v2.3.1.pdf) |
| Public source workbook (49 worksheets) | [Additional file 1](supplementary/Additional_file_1_Source_Data.xlsx) |
| Supplementary methods | [Word](supplementary/Additional_file_2_Methods.docx) · [PDF](supplementary/Additional_file_2_Methods.pdf) |
| Supplementary figures and legends | [Word](supplementary/Additional_file_3_Figures.docx) · [PDF](supplementary/Additional_file_3_Figures.pdf) |
| Main figures and Supplementary Figures 1 and 2 | `analysis_update/figures/` |
| Supplementary Figures 3 and 4 | `strengthening_v41/figures/` |
| Revised Table 3 | `tables/Table_3_v2.3.1.docx` and `.csv` |
| Aggregate plot inputs, including revised Figures 4A, 5A and 5B | `analysis_update/source_data/` |
| Original source access records | `metadata/original_input_files.csv` and the workbook |
| Changes and interpretation | [Release notes](metadata/RELEASE_NOTES_v2.3.1.md) |

The v2.3.1 files above supersede the older v41/v42 manuscript and supplementary files retained for provenance. The revised nerve interferon entries distinguish published atlas markers from the source polyneuropathy comparison. Figure 5A reports evidence types rather than a cross-tissue evidence ranking. CSF protein identification counts are not used as evidence of concentration increases.

## Reproduction

```bash
pip install -r analysis_update/requirements.txt
python analysis_update/analysis/run_v231_figures.py
python analysis_update/analysis/13_verify_public_summaries.py
```

The current plotting entry point recreates Figures 3–5 and standalone revised panels from the included aggregate tables. SVG text is editable; PDF fonts are embedded. Supplementary Figures 3 and 4 retain their aggregate plotting workflow in `strengthening_v41/`. Earlier figure entry points route revised panels to the current plotting functions. Full donor analyses require the original matrices and intermediate reconstruction documented in Additional file 2.

## Interpretation and distribution

Blood pooled effects remain non-significant. The macrophage Fc association weakens under correction across all 80 nerve comparisons (q = 0.087) and composition adjustment. No independent external validation or new functional experiment is claimed.

Only aggregate results and source access documentation are redistributed. Donor expression and clinical intermediates, the private author archive and manuscript Additional file 4 are excluded from this public repository. No original acquisition dates have been invented. The manuscript is a revised research draft; final journal and institutional requirements remain the authors' responsibility.

The concept DOI links the complete version series. The GitHub release page identifies the corresponding version-specific Zenodo record when available. Code is MIT licensed. Project-derived data and figures are CC BY 4.0 under `LICENSE-DATA`. External materials retain their original terms.
