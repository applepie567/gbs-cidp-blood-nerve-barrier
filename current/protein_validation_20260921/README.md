# PXD056286 protein quantification and data usability

This extension provides the completed extraction of stored nerve protein intensities, their instrument run mapping and a data usability assessment. The inspected public sources did not establish clinical specimen or diagnosis links, and the saved protein filter retained only FCGRT from the original 11-gene Fc module. These data could not support an independent disease comparison of the module. The quantitative exports and supporting records are retained for reuse.

The extraction and coverage audit was completed on 21 September 2026 after the manuscript findings were known. The Fc gene list comes from the preceding manuscript analysis; it was not selected to fit this protein dataset. This work is exploratory, with no preregistration and no disease hypothesis tests.

## Exported data

| File in results | Content |
| --- | --- |
| PXD056286_run_mapping.csv | All 19 instrument runs, native study and file IDs, quantitative columns, stored normalization factors, filename codes and six unavailable clinical metadata fields |
| PXD056286_raw_abundances.csv | All 2,356 protein records with stored quantitative values, before normalization, across F1 to F19 |
| PXD056286_normalized_abundances.csv | The corresponding 2,356 records with the normalization already stored in the source database |
| PXD056286_raw_stored_filter.csv and PXD056286_normalized_stored_filter.csv | The 1,460 quantified records retained by the saved protein filter |
| PXD056286_protein_QC.csv | Original identifiers, explicit FASTA gene annotations, master status, identification confidence, peptide support, group IDs and retained-filter flags |
| PXD056286_full_inventory.csv | All 18,630 protein records, including records without quantitative values |
| PXD056286_Fc_coverage.csv | Coverage and identification constraints for each of the original 11 Fc genes |
| normalization_verification.csv and validation_summary.json | Verification of all 34,107 available raw/normalized pairs and the overall extraction counts |

The manuscript source workbook (`../supplementary/Additional_file_1_Source_Data_strengthened.xlsx`) contains the run mapping, full quantitative matrices, protein QC and Fc coverage. The complete inventory, filtered matrices and source evidence are supplied here to keep the workbook manageable. Protein identifiers are strings, including signed 64-bit UniqueSequenceID values. Preserve them as text in Excel and when using JavaScript. They are protein record keys, not participant identifiers.

Values retain the source abundance scale. A stored missing-value flag produces an empty CSV or Excel cell. No missing value was set to zero, imputed or log transformed. No technical runs were averaged and no participant groups were inferred. The normalized matrix uses the normalization already stored by Proteome Discoverer; we did not estimate new normalization coefficients.

## What has been verified

Dynamic-column metadata link Abundances and AbundancesNormalized to distribution maps 49 and 53. Both maps contain the ordered labels F1: Sample to F19: Sample. StudyInformation and WorkflowInputFiles independently link these file IDs to the raw instrument filenames. For every available raw/normalized pair, the normalized value matches the raw value multiplied by the run coefficient stored in ProcessingNodeCustomData. The maximum relative discrepancy is 6.69e-16 or less.

The saved pdResultView filter requires Master Protein (IsMasterProtein 0), ExcludedBy -1, High protein FDR confidence (code 3), and GroupUniquePeptidesCount at least 2. This is the source's saved filter, not a newly chosen target-specific threshold. Its reconstructed ID set exactly matches all 1,470 entries in TargetProteins_Subset_4; 1,460 have quantitative values. GroupUniquePeptidesCount and UniquePeptidesCount are distinct fields. The publication reports 1,411 proteins, so unrecorded downstream selection or aggregation remains unresolved. We do not claim to have reproduced that final published matrix.

FCGRT, FCGR3A, FCGR3B, LYN and HCK have stored quantitative values. Only FCGRT passes the saved filter and it is observed in 12 runs. FCGR3A and FCGR3B each have only one group-unique peptide. LYN and HCK are low-confidence master candidates and share protein group 284. The gene list is FCGR1A, FCGR2A, FCGR2B, FCGR3A, FCGR3B, FCGRT, FCER1G, TYROBP, SYK, LYN and HCK. Absence of a quantitative record is not evidence of biological absence. We used the original GeneSymbol field and explicit GN annotations in the stored FASTA description to audit coverage; no new sequence search or ambiguous peptide reassignment was performed. The raw GeneSymbol field remains unchanged, including its missing value for FCGR3B accession O75015.

## Clinical mapping remains unresolved

The publication describes nine CIDP cases and two non-diseased controls in its proteomic comparison. The deposited analysis has 19 instrument runs and 13 distinct filename codes. N241, VNP132 and VNP141 each occur three times. Neither the 19 runs nor the 13 filename codes can be treated as a confirmed biological sample count.

StudyInformation SampleGroup and ReplicateGroup fields are empty for all runs, and the AnalysisDefinition contains no study factors or sample factor values. The inspected pdResultView and msfView metadata do not supply a clinical crosswalk. The first 65,536 bytes of each of the 19 raw files contain generic AH_CIDP_Nerv labels; these are not verified diagnoses or Figure 3 Patient 1 to 9 labels. The associated Figure 3 and supplement captions do not link those patient labels to raw filenames. This audit covers the specified metadata and raw file headers, not every byte of the raw acquisition files.

Columns ending in _author in the run mapping remain intentionally blank because the inspected public sources did not establish the specimen ID, diagnosis, technical replicate group, paper sample label, inclusion in the published analysis or a mapping evidence reference. The suffix is retained from the original export schema. Repeated filename codes and clustering of protein abundances cannot establish these links.

The retained protein coverage is insufficient to reproduce the original 11-gene Fc module. FCGRT alone is a different endpoint. Whole nerve protein abundance also does not isolate macrophages and cannot resolve cell-composition confounding of the RNA result. Participant independence from the discovery atlas was not established. No disease contrast or disease P value was generated from these runs. This completed usability assessment is not a positive or negative replication result.

## Reproduction

Use Python 3 with pandas and NumPy. The original execution used Python 3.12.14, pandas 2.2.3 and NumPy 2.3.5. The retrieval and binary decoding scripts otherwise use the Python standard library.

Run from this directory, sequentially and wait for each command to finish:

```bash
python extract_mapping_metadata.py
python extract_quantification.py
python inspect_raw_headers.py
python assess_exports.py
```

Binary cached byte ranges and raw-header fragments are omitted from this compact update. Run the extraction commands with internet access to retrieve the required ranges; downloading the complete 16,144,965,632-byte result database is not required. The aggregate exports and source metadata snapshots are included. Source metadata JSON files are snapshots; extract_mapping_metadata.py preserves existing snapshots. To independently refresh them, work in a separate copy and remove the relevant snapshot and cache files before running the script with internet access. All output tables are regenerated by the extraction and assessment scripts. The source_run_mapping.csv was taken from the preceding metadata audit and is checked against both newly retrieved native study and workflow tables on every assessment run.

Abundance BLOBs are decoded as repeated little-endian double/presence-byte pairs; present flags are checked and missing flags retained. The lower-level SQLite table reader was compared with Python sqlite3 on multi-page tables including signed values, UTF-8 text and overflow BLOBs. The full extraction also checks unique protein keys, vector lengths, shared missingness, exact native filter membership, native run links and stored normalization identities. The public source is read without changes.

access/*_ranges.json records byte offsets, sizes and SHA256 hashes of the corresponding cached database segments. raw_header_review.json records the source URL, byte count and SHA256 of each raw header fragment. The repository-root MANIFEST_SHA256.json records all files delivered in this update. It does not claim a hash of the entire remote result database.

## Sources

- Stascheit and colleagues, associated nerve study: https://doi.org/10.1007/s00401-025-02936-w
- PRIDE project: https://www.ebi.ac.uk/pride/archive/projects/PXD056286
- ProteomeXchange record: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD056286
- Public source directory: https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD056286/
- Main result file: QExactiveHF02_18988_N1530.pdResult; saved views: QExactiveHF02_18988_N1530.pdResultView and QExactiveHF02_18988_N1530.msfView
- SQLite file format used for the range reader: https://www.sqlite.org/fileformat.html

These files are included in the prepared v2.4.0 manuscript update. No repository push, release publication or new archive deposit was performed while preparing this package.
