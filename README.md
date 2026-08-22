# Myeloid–blood–nerve barrier–Schwann cell circuit in GBS and CIDP

Reproducibility repository for the manuscript:

> **Cross-platform multi-omics analysis maps a myeloid–blood–nerve barrier–Schwann cell circuit across Guillain–Barré syndrome and CIDP**

This project integrates acute Guillain–Barré syndrome (GBS) blood transcriptomics, paired longitudinal serum proteomics, chronic inflammatory demyelinating polyradiculoneuropathy (CIDP) sural-nerve single-nucleus transcriptomics, and a normal human blood–nerve barrier (BNB) reference.

## Anatomical terminology

The BNB is the endoneurial microvascular endothelial barrier of peripheral nerve. It is not referred to as the blood–brain barrier (BBB). The perineurium is a distinct, nonvascular diffusion barrier and is treated separately.

## Data resources

Seven human resources are analyzed directly:

1. **GSE211225** — whole-blood transcriptomics: 6 acute GBS, 10 post-acute GBS, and 6 controls.
2. **GSE31014** — peripheral-blood leukocyte microarray: 7 GBS and 7 controls.
3. **GSE304871** — early untreated AIDP sorted CD11b+, CD4+, and CD8+ RNA-seq.
4. **PRJNA1293757** — PBMC single-cell RNA-seq: 3 treatment-naive AIDP and 2 controls.
5. **GBS-Proteomics** — paired acute-to-one-year SomaScan measurements from 20 patients and 15 controls.
6. **GSE285983** — 365,708 nuclei from 37 sural nerves, including 9 CIDP and 11 CIAP donors.
7. **GSE107574** — cultured human endoneurial endothelial cells and laser-captured endoneurial microvessels used as a normal BNB reference.

OEP002315/OEP002701 provides author-reported monocyte pathway support. GSE133750 rat experimental autoimmune neuritis is used as an animal mechanism tier. Accession links and evidence roles are listed in `metadata/data_source_urls.csv`.

## Repository contents

- `analysis/` — numbered analysis workflows.
- `config/` — dataset and contrast specifications.
- `scripts/` — public-data retrieval helpers.
- `data/README.md` — data provenance, expected paths, and redistribution boundaries.
- `results/` — machine-readable derived results.
- `source_data/` — panel-level source data for all main figures and tables.
- `figures/` — publication figures in PNG and vector PDF formats.
- `docs/study_protocol.md` — extended computational protocol.
- `metadata/` — dataset registry and release metadata.

Raw high-dimensional matrices and participant-level clinical records are not redistributed. Obtain primary data from the originating repositories under their respective terms.

## Reproduction

Create an isolated Python environment and install the recorded dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run the numbered workflows in order after acquiring the source datasets and setting the expected paths described in `data/README.md`:

```bash
python analysis/01_targeted_pseudobulk.py
python analysis/02_gbs_sorted_bulk.py
python analysis/04_prjna1293757_pseudobulk.py
python analysis/05_gbs_bulk_crossplatform.py
python analysis/06_gbs_longitudinal_proteomics.py
python analysis/07_bnb_ean_external_validation.py
python analysis/08_integrate_upgrade_figures.py
```

The participant or tissue donor is the unit of inference. Cells and nuclei are not treated as independent biological replicates.

## Citation

Please cite the archived release using the metadata in `CITATION.cff`. The permanent Zenodo DOI will be inserted after the release metadata and author details have been approved.

## Licenses

- Analysis and plotting code: MIT License (`LICENSE`).
- Repository-authored aggregate source data and documentation: CC BY 4.0 (`LICENSE-DATA`).
- Primary datasets remain subject to the licenses and access conditions of their source repositories.

## Contact

Corresponding-author email addresses and institutional details will be added after author confirmation and before the public release.
