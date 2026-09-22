# Reproducing the version 41 additions

All new tests are exploratory and follow an already observed association.

`strengthen_fc.py --source ORIGINAL_RESULTS --out RESULTS` reads the original donor expression table, module table and GSE285983 pseudobulk JSON. The original results directory contains a `tables` subdirectory. These donor tables can be reconstructed with the archived `01_cidp_nerve_pseudobulk.py` and all 37 original matrices. They are not included in the public aggregate package.

`download_inputs_v41.py --manifest results/new_input_manifest.csv --out LOCAL_INPUTS` retrieves the 20 matrices used in the new annotation analysis, both annotation files and all eight spatial matrices, with SHA256 validation. This downloader does not acquire the other 17 nerve matrices needed to rebuild the original 37 donor standardisation.

`subtype_analysis.py --root LOCAL_INPUTS --original-results ORIGINAL_RESULTS` expects the output from strengthen_fc.py at LOCAL_INPUTS/results. It reaggregates the downloaded original matrices. The author archive includes a validated targeted pseudobulk checkpoint, which supports the additional option `--reuse-checkpoint` without redownloading the matrices.

`plot_strengthening.py --results results --out regenerated_figures` regenerates Supplementary Figures 3 and 4 from the aggregate tables. Times New Roman must be installed. Adjust the two font installation paths in this plotting script to match a different machine.

Python 3.12.13, numpy 2.3.5, pandas 2.2.3 and scipy 1.17.0 were used. h5py reads original counts. matplotlib produces the figures. Exact installed versions are also recorded in environment_versions.json.

The earlier analysis_update folder retains the version 39 CSF correction and earlier figure code. Its cross compartment evidence table has been updated to include the version 41 Fc caveat. The existing Zenodo v2.1.0 DOI is a preceding release and does not contain these new files.

Public donor omission results are ranges rather than individual donor rows. The private author archive retains the complete omission table and other donor intermediates. External original materials remain subject to their source terms.
