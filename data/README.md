# Data retrieval and provenance

Raw high-dimensional matrices are intentionally excluded from the submission
archive. Derived target-gene tables, sample maps and statistical outputs are
included under `results/tables/`.

## Direct human resources

- **GSE211225**: processed whole-blood matrix and metadata from GEO. Acute and
  post-acute participants are different groups.
- **GSE31014**: processed GPL96 leukocyte matrix from GEO.
- **GSE304871**: normalized sorted-cell matrices from GEO.
- **PRJNA1293757**: processed PBMC single-cell objects. Use
  `scripts/download_prjna1293757_processed.sh`; keep sample identity and aggregate
  module scores to the participant level.
- **GBS-Proteomics**: public SomaScan data and scripts from
  <https://github.com/NeRveLabBCN/GBS-Proteomics>. Match acute and one-year samples
  using the repository matching identifier. Label all results as preprint evidence.
- **GSE285983**: download `GSE285983_metadata_all.csv.gz` and `GSE285983_RAW.tar`
  from GEO and extract H5 matrices into `data/raw/GSE285983/h5/`. Aggregate counts
  by donor and compartment before inference.
- **GSE107574**: processed FPKM matrix from GEO. Cultured endothelial preparations
  and laser-captured microvessels are descriptively summarized and must not be
  tested as disease groups.

## External and animal resources

- **OEP002315/OEP002701**: only author-supplied monocyte GSEA tables from the
  accessible supplementary package are used. Do not describe these as a direct
  individual-level reanalysis.
- **GSE133750**: rat sciatic-nerve RNA-seq at control, early, peak and late EAN
  (n=3/group). Human targets are mapped to rat orthologues; Cxcl1/Cxcl2 with
  Cxcr1/Cxcr2 represent the rat CXCL8-axis surrogate.

Official accession and publication URLs are recorded in
`results/tables/data_source_urls.csv` and Additional file 1.

## Anatomical terminology

BNB denotes peripheral nerve endoneurial microvascular endothelium and associated
pericytes. BBB is reserved for the central nervous system. The perineurium is
analyzed separately.
