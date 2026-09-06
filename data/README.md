# Original inputs and local intermediates

Original filenames, accession numbers, source URLs and processing scripts are listed in `metadata/original_input_files.csv` and the `Input_files` sheet of Additional file 1. Sources are also summarized in `config/datasets.yaml` and `metadata/data_source_urls.csv`.

Primary pipelines expect inputs under `data/raw/GSE211225/`, `data/raw/GSE31014/`, `data/raw/PRJNA1293757/` and `data/raw/GSE285983/`. Follow the exact filenames and archive-member paths in the input register. GSE107574 is a published BNB reference; a separately processed file is not documented in the supplied package.

Participant and donor expression/clinical intermediates are regenerated locally and omitted from the public package. See `metadata/data_reconstruction.csv` for omitted files and reconstruction requirements. The supplied genetic scripts start from a precomputed donor-cell table and do not provide its raw-matrix reconstruction.

The public workflow uses aggregate tables and published coefficients. Primary pipelines were not rerun for this documentary update. Original download dates and input checksums were not recorded, so none are invented here.
