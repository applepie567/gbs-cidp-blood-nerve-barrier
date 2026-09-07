# GBS and CIDP immune compartmentalization

**Immune compartmentalization across blood, cerebrospinal fluid and peripheral nerve in Guillain-Barré syndrome and CIDP**

Version **2.1.0** accompanies manuscript **v37**, with final figures, aggregate results and original-file access information.

## Current materials

| Material | Location |
|---|---|
| Manuscript v37 | `docs/GBS_CIDP_compartmentalization_v37_data_availability_2026-09-06.docx` |
| Five main figures and two supplementary figures | `figures/` |
| Three manuscript tables | `tables/` |
| Additional file 1: 30-sheet public workbook | `source_data/Additional_file_1_source_data_public.xlsx` |
| Additional file 2: reproducibility and source appendix | `docs/Additional_file_2_reproducibility_and_sources_public.docx` |
| Nineteen original-input/source records | `metadata/original_input_files.csv` |
| Local reconstruction requirements | `metadata/data_reconstruction.csv` |
| Figure/table source index | `source_data/Figure_table_index.csv` |
| Recorded reproduction environment | `metadata/REPRODUCTION_ENVIRONMENT.json` |

## Run with the public aggregate tables

Use Python 3.12 and the recorded dependencies:

```bash
python -m pip install -r requirements-reproduction.txt
python analysis/run_v2_release.py
python analysis/render_current_figures.py
python tests/validate_release.py
```

The default workflow recalculates blood random-effects synthesis and leave-one-cohort-out estimates, CSF coefficient pooling, the interval for the published 52/55 biopsy count, and multiple-testing corrections from supplied aggregate P values. It compares these calculations with the frozen results. It checks all workbook sheets against their CSV sources and all seven figure PNGs against manuscript v37. The five main figures can be rendered from the supplied aggregate inputs.

Genetic localization summaries and bootstrap intervals, nerve contrasts and correlation estimates remain frozen results. Checking their P-value corrections does not rerun donor-level tests, correlations or bootstrap resampling. Figure 4A uses its verified raster interior in `analysis/assets/figure4_localization_reference.png`; its labels and remaining panels are drawn by code. Supplementary PNGs are preserved from the manuscript with their aggregate numerical inputs. Other rendering environments may produce pixel differences.

## Recalculate from local donor tables

Individual donor/participant intermediates and JSON files containing sample-level measurements are omitted from this version. The workbook retains 28 original aggregate/documentation sheets and adds `Input_files` and `Data_reconstruction`; it excludes `CIDP_donor_heterogeneity` and `Genetic_donor_celltype`.

Download original inputs using the source records. Primary pipelines are `01_cidp_nerve_pseudobulk.py`, `02_prjna1293757_pseudobulk.py` and `03_gbs_blood_crosscohort.py`; broader dependencies are in `requirements.txt`. GSE285983 processing reconstructs targeted expression, module scores and cell fractions locally. The five selected CIDP module scores must then be pivoted into the donor-by-module correlation input.

The supplied genetic compilation and bootstrap scripts require the precomputed `Genetic_donor_celltype.csv`. They do not reconstruct that table from raw matrices. Exact end-to-end genetic reproduction requires the original intermediate or a separately implemented and validated reconstruction. After all required local donor inputs exist, run:

```bash
python analysis/run_v2_release.py --with-local-donor-tables
```

This option checks required files before beginning. Local individual-level outputs are excluded by `.gitignore`. Primary matrices were not rerun for this update. Original download dates and input checksums were not recorded, so no retrospective values are assigned.

## Evidence and interpretation

Blood synthesis uses three cohorts and modified Hartung–Knapp scaling bounded below by one. CIDP nerve analyses use donors as the units of inference; reported module contrasts are unadjusted for age, sex and center. CSF pooling uses published adjusted coefficients rather than new participant-level CSF observations. Genetic association estimates are attributed to Du et al. 2024, and tissue observations to Stascheit et al. 2025. Cross-compartment comparisons use independent cohorts and do not establish temporal progression from GBS to CIDP or causal drug targets.

## Versions and citation

The current published archive is [v2.1.0, DOI 10.5281/zenodo.22561445](https://doi.org/10.5281/zenodo.22561445), published on 7 September 2026. The [GitHub v2.1.0 release](https://github.com/applepie567/gbs-cidp-blood-nerve-barrier/releases/tag/v2.1.0) provides the identical public ZIP. The archived package and tag preserve commit `9682d3374020af5cc6c236613abf3355630024ef`; this subsequent citation update records the issued DOI.

The preceding archive, [v2.0.0, DOI 10.5281/zenodo.22226674](https://doi.org/10.5281/zenodo.22226674), predates the extended analyses and revised attachments.

Older tags and archives remain historical records. This update removes obsolete preprint outputs and individual-level copies from the current branch without rewriting repository history or earlier archives. See `metadata/RELEASE_NOTES_v2.1.0.md` and `metadata/PUBLIC_DISTRIBUTION.json`.

Code is MIT licensed. Project-derived data and figures are CC BY 4.0 under `LICENSE-DATA`. External datasets and publications retain their original terms and citation requirements.
