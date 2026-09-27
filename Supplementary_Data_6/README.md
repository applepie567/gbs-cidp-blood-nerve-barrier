# Nerve macrophage marker analysis

This supplement extends the GSE213455 whole nerve analysis. It does not provide independent cell-resolved CIDP versus CIAP replication.

## Result

The reference contains 365,708 author-retained nuclei from 37 GSE285983 donors, aggregated into 13 broad lineages. The 20-marker index excludes all 11 Fc outcome genes. Leave-one-donor-out Spearman correlations are 0.933381 with the observed macrophage nucleus fraction and 0.950450 with the macrophage share of captured counts. Leaving out entire centres gives 0.911332 and 0.927217, respectively. All within-centre checks remain positive and the recorded feasibility criteria are met.

In the independent four-CIDP versus nine-vasculitic-neuropathy cohort, the fixed marker index correlates with the Fc score at rho 0.906593. The original difference is -0.836553 (HC3 95% CI -1.763592 to 0.090486). Adding the index gives 0.208146 (-0.915957 to 1.332250, model P 0.688623). Mean-probe and donor-omission sensitivities are supplied in full. The conditional result concerns marker expression and does not separate cell abundance from cell state. The index is not an absolute cell fraction, and same-platform atlas calibration does not validate cross-platform accuracy.

## Reproduce the statistical analysis

From this directory, run:

```
python scripts/analyse_composition.py
python scripts/make_figure.py
Rscript scripts/verify_models.R .
```

The analysis uses the included all-gene donor-lineage counts, metadata, external series matrix, GPL13369 annotation, and unchanged external Fc scores. It needs no large matrix download. Gene ranks use a fixed common background of 16,633 genes, excluding the 11 outcome genes. Marker selection is rebuilt within every held-out training set. `results/all_marker_selections.csv` retains every selected gene and selection statistic. The NPZ array has rows matching `atlas_donor_lineage_metadata.csv`, and its `genes` entry identifies columns. Counts include all deposited gene symbols after summing duplicates. The original 11-gene Fc reconstruction is exact for every donor.

To repeat aggregation from the public source matrices, install h5py or use an existing compatible HDF5 library and run:

```
python scripts/download_reference.py
python scripts/build_reference.py
python scripts/analyse_composition.py
```

The download script retrieves and verifies the 37 matrices and full original nucleus metadata. Source URLs and hashes are in `source_manifest.csv`. The compact `atlas_nucleus_annotation.csv.gz` records the retained nucleus labels used here. The full original metadata and large HDF5 files are not redistributed in this package.

## Figure and table sources

- S12A: `results/donor_holdout_predictions.csv`
- S12B: `results/centre_holdout_predictions.csv`
- S12C: `results/external_donor_scores.csv`
- S12D and Table S6: `results/external_conditional_models.csv`
- Additional calibration metrics: `results/feasibility_metrics.csv`
- Marker identities: `results/frozen_macrophage_markers.csv`
- Omission sensitivity: `results/external_leave_one_out.csv`

All points are donors, not individual nuclei. The original Fc score is unchanged. Conditional models use the existing cohort score scale, including during omission checks. Model intervals and P values are separate from the original exact rank P value.

## Provenance and availability

Heming et al., Nature Communications 2025, doi:10.1038/s41467-025-62964-8, GSE285983.
Cerri et al., Frontiers in Oncology 2022, doi:10.3389/fonc.2022.974751, GSE213455 and GPL13369.

This is an exploratory extension specified after the original external disease result was known. Its marker and feasibility rules were written before deriving the new marker results. Public data retain their source attribution and terms. The previous public Zenodo v2.4.0 archive does not contain this analysis. No new DOI is claimed.
