# Immune compartmentalization across blood, CSF and peripheral nerve in GBS and CIDP

This repository contains the Python analysis workflow, derived source data, figures and reproducibility materials for:

**Immune compartmentalization across blood, cerebrospinal fluid and peripheral nerve in Guillain-Barré syndrome and CIDP**

Version 2.0.0 restructures the study around three biological compartments:

- acute GBS blood transcriptomes
- published GBS cerebrospinal-fluid proteomic evidence
- CIDP sural-nerve single-nucleus profiles and blood-nerve barrier localization

Published CIDP genetic evidence is anchored to donor-resolved peripheral-nerve and blood-nerve barrier cell types. The genetic analysis prioritizes gene-cell combinations for functional validation and is not presented as direct drug-target discovery.

## Repository contents

- `analysis/01_cidp_nerve_pseudobulk.py` performs donor-level targeted pseudobulk analysis of GSE285983.
- `analysis/02_prjna1293757_pseudobulk.py` performs biological-sample pseudobulk analysis of PRJNA1293757.
- `analysis/03_gbs_blood_crosscohort.py` analyzes GSE211225 and GSE31014.
- `analysis/04_compile_published_csf_evidence.py` compiles the published CSF evidence while preserving its original inference level.
- `analysis/05_compile_genetic_cell_localization.py` exports CIDP genetic evidence and donor-resolved cell localization.
- `analysis/06_build_cross_compartment_tables.py` builds the blood-CSF-peripheral-nerve evidence map.
- `analysis/07_export_workbook_sheets.py` exports every source-workbook sheet as CSV.
- `analysis/run_v2_release.py` runs the release-level compilation and validation steps.
- `source_data/` contains the complete machine-readable source workbook and individual CSV sheets.
- `figures/` contains the five publication figures at their original embedded resolution.
- `tables/` contains machine-readable versions of Tables 1–3.
- `docs/` contains the manuscript and reproducibility appendix.
- `metadata/` contains dataset URLs, release notes, checksums and validation results.

## Reproduce the release-level outputs

Install Python 3.12 or later and the listed dependencies:

```bash
python -m pip install -r requirements.txt
python analysis/run_v2_release.py
```

This command exports the published CSF evidence, genetic localization tables, cross-compartment tables and all workbook sheets, then validates the release contents.

The raw-data scripts require public repository files under `data/raw/`. Download instructions and expected paths are provided in `data/README.md`. Large primary sequencing and mass-spectrometry files are not redistributed.

## Statistical units

Blood cohorts are analyzed at the participant or biological-sample level. CIDP nerve analyses aggregate nuclei to donor-level pseudobulk profiles. CSF evidence retains the inference level reported in each publication and does not convert group-level findings into participant-level observations.

## Data and code availability

The current repository is:

https://github.com/applepie567/gbs-cidp-blood-nerve-barrier

The historical v1.0.0 archive is available at:

https://doi.org/10.5281/zenodo.22067879

The published v2.0.0 archive is available at:

https://doi.org/10.5281/zenodo.22226674

## Licensing

Code is released under the MIT License. Derived and aggregated source data are released under CC BY 4.0. The original public datasets remain subject to their repository and study-specific terms.
