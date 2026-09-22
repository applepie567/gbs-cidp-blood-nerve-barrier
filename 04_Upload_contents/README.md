# Immune recruitment signatures and cellular context in inflammatory neuropathies

Prepared manuscript update **v2.4.0 — 22 September 2026**.

This package is a self-contained update for [the existing repository](https://github.com/applepie567/gbs-cidp-blood-nerve-barrier). The current manuscript and reproducibility materials are in **`current/`**. Existing versioned directories in the repository can remain as historical material. This archive has not been pushed to GitHub or deposited on Zenodo.

## Start here

- [Current manuscript](current/manuscript/Manuscript_strengthened.docx)
- [Figure source workbook](current/figure_source_data/Figure_source_data_v2.4.0.xlsx)
- [Per-panel CSV index](current/figure_source_data/INDEX.csv)
- [Current figures](current/figures/)
- [Current figure legends](current/metadata/current_figure_legends.csv)
- [Public proteomics exports and feasibility assessment](current/protein_validation_20260921/README.md)
- [Chinese update instructions](UPLOAD_GITHUB_zh.md)

## Reproduce all figures from supplied source tables

```bash
python -m pip install -r requirements.txt
python current/analysis/reproduce_figures.py --out rebuilt_figures
python current/analysis/verify_release.py
```

The supplied PNGs are identical to the images embedded in the current Word files. Figure 1, Figure 3 and Figure 5 incorporate the final removal of redundant in-figure sentences. Relevant statistical detail is retained in the legends and source data. Figure S2 also has an editable reconstruction; its layout can differ from the retained author-checked PNG, while the published count and observations are unchanged. Font substitution can alter layout; use Times New Roman for the closest match, or the documented serif fallback. Fonts are not redistributed.

The CSVs contain the values underlying every plotted point, interval, heatmap and published count. They include anonymous donor-level coordinates where used in Figure 5, aggregate expression and effect estimates elsewhere, and published CSF/genetic coefficients. Figure 1 is a schematic. These tables are **figure source data**, not a substitute for original sequence reads or participant-level CSF measurements.

## Analysis and provenance

`current/analysis/` contains the current exploratory extension and plotting entry point. `current/legacy_reference/` supplies the aggregate inputs for retained analyses; `current/original_analysis/` preserves the earlier pipeline. `current/followup_20260921/` contains the graft-reference and 91-test follow-up. Raw-input URLs and hashes are in `current/source_data/public_input_manifest.csv`; see [reanalysis instructions](current/REPRODUCIBILITY.md).

Distinct testing families remain separate: within-cell module correction, expanded nerve families, the nine new programme comparisons, and the 17-test clinical sensitivity. The newly prepared tables expose the current 91-test correction while retaining historical 80/89-test records for provenance.

Public PXD056286 quantification exports are available. Reliable run-to-donor-to-diagnosis links and the final disease-labelled analysis matrix remain unavailable, so independent nerve validation was **not performed**. The feasibility assessment is not a negative validation result. Remote database byte caches are retrieved on demand and omitted from this compact package.

## Versions, archive and citation

The preceding live release was `V2.3.1` when this package was prepared. The manuscript's existing version-specific archive citation is retained in the frozen Word file. After publishing this update, replace that citation with the actual new version and DOI. The known series/concept record is [10.5281/zenodo.22803504](https://doi.org/10.5281/zenodo.22803504); no new version DOI has been assigned here.

Code: MIT. Project-derived data and figures: CC BY 4.0. External materials retain their original licences. See `LICENSE`, `LICENSE-DATA`, `CITATION.cff`, and `current/metadata/release.json`.
