# Immune gene and protein profiles across blood, cerebrospinal fluid and peripheral nerve in inflammatory neuropathies

Release **v2.3.0** accompanies the Human Genomics submission manuscript **v42** (7 September 2026).

Archive: [DOI 10.5281/zenodo.22641099](https://doi.org/10.5281/zenodo.22641099). GitHub: [v2.3.0 release](https://github.com/applepie567/gbs-cidp-blood-nerve-barrier/releases/tag/v2.3.0).

| Material | Location |
|---|---|
| Current manuscript with five embedded figures and three editable tables | [Word](GBS_CIDP_Human_Genomics_v42.docx) and [PDF](GBS_CIDP_Human_Genomics_v42.pdf) |
| Current graphical abstract | [PNG](graphical_abstract/GBS_CIDP_Human_Genomics_graphical_abstract.png) |
| Historical v41 text manuscript | [Word](GBS_CIDP_v41_strengthened.docx) and [PDF](GBS_CIDP_v41_strengthened.pdf) |
| Five main figures and Supplementary Figures 1 and 2 | `analysis_update/figures/` |
| Supplementary Figures 3 and 4 | `strengthening_v41/figures/` |
| Three tables | `tables/` |
| Public workbook, 49 worksheets | `supplementary/Supplementary_Data_1_v41.xlsx` |
| Supplementary methods and figure legends | `supplementary/` |
| Original aggregate figure source tables | `analysis_update/source_data/` and `analysis_update/results/tables/` |
| Added aggregate source tables and new input manifest | `strengthening_v41/results/` |
| Original input file names, databases and source URLs | `metadata/original_input_files.csv` and the workbook |
| Changes and interpretation | `metadata/RELEASE_NOTES_v2.3.0.md` |

The v42 manuscript supersedes the v41 manuscript files retained in this repository. The v41 supplementary workbook and figures continue to accompany v42. Supplementary Methods retains its v41 filename, with the archive links updated. Analysis results are unchanged from v2.2.0.

## Reproduction

The scripts and frozen aggregate results are included. See [the v41 instructions](strengthening_v41/README.md) and Supplementary Methods for required original inputs and limitations. Supplementary Figures 3 and 4 can be regenerated from the aggregate tables:

```bash
python strengthening_v41/plot_strengthening.py --results strengthening_v41/results --out regenerated_figures
```

The plotting script requires Times New Roman and its recorded Python dependencies. Adapt font paths for a different machine. Earlier figure scripts are under `analysis_update/analysis/`, including `run_v39_figures.py`. The legacy `run_v2_release.py` refers to the older release layout and is retained for provenance, not as the current package entry point.

Full donor analyses require local reconstruction from original data. The new downloader retrieves the 20 matrices used for annotation reaggregation, the annotation files and eight spatial matrices. The other 17 nerve matrices are also needed to reconstruct the original 37 donor standardisation. Full genetic localization reconstruction requires the additional intermediate described in Supplementary Methods.

## Interpretation and distribution

The added tests are exploratory. The Fc association weakens with broader multiplicity correction and composition adjustment. No independent external cohort validation or functional experiment was completed. The spatial panel does not cover the Fc genes.

Only aggregate results and source access documentation are redistributed. Donor expression and clinical intermediates and the private author archive are excluded. Earlier source acquisition dates were not recorded and have not been invented.

Prior release records include [the manual v2.1.0 archive](https://doi.org/10.5281/zenodo.22561445) and [its GitHub mirror](https://doi.org/10.5281/zenodo.22563142). The concept DOI [10.5281/zenodo.22226672](https://doi.org/10.5281/zenodo.22226672) represents all versions. Cite the specific version used for an analysis.

Code is MIT licensed. Project derived data and figures are CC BY 4.0 under `LICENSE-DATA`. External materials retain their original terms.
