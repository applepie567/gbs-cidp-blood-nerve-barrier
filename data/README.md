# Primary-data locations

Large primary files are not stored in this repository. The analysis scripts expect the following local structure:

```text
data/raw/
├── GSE211225/
├── GSE31014/
├── PRJNA1293757/
├── GSE285983/
└── GSE107574/
```

The source repository URLs and analytical roles are defined in `config/datasets.yaml` and `metadata/data_source_urls.csv`.

PXD002911 and the additional published CSF studies are compiled at their reported inference level in the source workbook. The release does not generate patient-level CSF observations from published group summaries.

