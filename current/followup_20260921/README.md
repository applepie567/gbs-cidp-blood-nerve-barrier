# Focused Fc receptor follow up

This exploratory follow up was completed on 21 September 2026 after the earlier
manuscript findings were known. It adds focused graft reference comparisons,
checks whether direct macrophage composition standardisation is supported, and
audits the public PXD056286 result database. It does not add an independently
validated disease effect. The existing public release was not changed.

## Findings

The original 11 gene macrophage score was reproduced for all 37 atlas donors to
within 1e-12. The prior CIDP minus CIAP difference remains -0.673. Its original
Mann Whitney P value, 0.0010874, is retained in the manuscript and expanded
correction. The 91 test sensitivity q value is 0.09896.

Relative to four graft reference donors, the mean score differences were -0.401
for CIDP and 0.271 for CIAP. The respective bootstrap 95% intervals were -0.901
to 0.050 and -0.171 to 0.652. Exact two sided Mann Whitney P values were 0.33007
and 0.34286. Both q values were 0.34286 across the two focused reference tests
and 1.0 in the 91 test sensitivity family. The small reference group, recruited
at one centre and comprising one male and three female donors, does not resolve
the direction of disease change relative to reference tissue. Nonsignificance
does not establish equivalence or a normal expression level.

The original pipeline had already calculated CIDP versus CTRL comparisons.
This follow up reconstructs that Fc result, adds CIAP versus CTRL and reports
both in a focused family. The 91 test sensitivity extends the manuscript's
previous 89 test family. It is not a correction across every analysis ever run
in the original pipeline. All tests remain exploratory. The exact P value for
the prior CIDP versus CIAP row in the audit table does not replace the original
asymptotic test or rescue its expanded significance.

None of 18 refined macrophage subtypes met the existing 20 nucleus threshold
in all 20 disease donors. Two CIDP donors and one CIAP donor had no eligible
fine subtype. Full cohort direct standardisation at this resolution would
require unobserved subtype expression to be estimated, so no standardised
effect was reported. This finding does not prove that all of the original
association is caused by composition. Existing covariate and subtype analyses
remain within-cohort sensitivity analyses.

PXD056286 contains protein abundance fields and 19 instrument runs. The
inspected StudyInformation table has no sample group or replicate group
assignments, and the AnalysisDefinition has no study factors or sample factor
values. The public article reports nine CIDP cases and two controls in its
proteomic subset. Run names alone do not establish how the deposited runs map
to those biological specimens. Protein quantities were therefore not decoded,
and no protein disease effects or coverage conclusions were estimated. A
verified specimen and diagnosis mapping plus a protein quantitative export
could support further tissue-level analysis. Donor independence and the
difference between whole tissue proteins and cell-specific RNA would still
need to be assessed.

## Reproduce the nerve follow up

Unpack the existing `GBS_CIDP_v41_author_archive_2026-09-07.zip`. Its
`03_private_reanalysis_workspace` contains the original donor pseudobulk
expression, archived module scores and refined subtype nucleus counts. Input
filenames, SHA256 values and software versions are listed in
`results/nerve_followup_summary.json`.

From this extension archive's root, replace `/path/to/archive` below with the
unpacked v41 archive root. The command regenerates aggregate results and writes
donor outputs separately. It does not download or publish any donor records.

```bash
python followup_20260921/nerve_fc_followup.py \
  --source /path/to/archive/03_private_reanalysis_workspace \
  --extension . \
  --out followup_20260921/results \
  --private-out local_followup_donor_outputs
```

The independent unit is the donor. The score retains standardisation across
all 37 atlas donors. Bootstrap intervals use 20000 resamples independently
within each comparison group, seed 20260921, with the score standardisation
held fixed. These intervals describe mean differences, whereas rank tests
compare score distributions. The inference table keeps both methods labelled.
No donor expression or clinical tables have been added to this attachment.

## Reproduce the protein metadata audit

Metadata snapshots and the exact cached byte ranges are supplied in
`protein_metadata`. The range manifest records source URLs, byte positions and
SHA256 values. The database itself was not modified. A metadata-only extraction
avoids downloading the approximately 16 GB result file.

```bash
python followup_20260921/read_sqlite_ranges.py --self-test
python followup_20260921/audit_proteome.py \
  --access followup_20260921/protein_metadata \
  --out followup_20260921/results
```

Add `--remote` to fetch metadata from the original PRIDE file. Supplied cached
ranges are reused. To independently fetch the ranges again, use a new empty
`--access` directory. The reader validates HTTP range responses and is checked
against Python sqlite3 using Unicode text, numeric values, blobs, interior
pages and overflow records. Only SQLite table b-trees are implemented. Protein
abundance BLOB decoding is outside this audit.

## Attachment changes

`results` contains complete aggregate reference contrasts, the 91 test family,
subtype support, group summaries and the protein feasibility audit. Blank q
values for the original contrast mean that it is outside the focused reference
family. Blank protein group fields reproduce missing source metadata.

The source workbook adds eight sheets and changes four README cells, the nerve
donor count and protein audit role in Table 1, and the Fc receptor cell in
Table 3. Earlier analysis sheets retain their numeric values.
`source_data/manuscript_tables.json` is updated to match the manuscript tables.
Existing figure source files are unchanged, because the plotted
estimates and significance symbols are unchanged. The manuscript captions
state the expanded correction scope. Detailed changes are listed in
`results/manuscript_followup_changes.json`. Authoring scripts are included for
provenance and use the workspace layout recorded in those scripts.

## Sources

- Nerve atlas https://doi.org/10.1038/s41467-025-62964-8
- Protein study https://doi.org/10.1007/s00401-025-02936-w
- Protein repository https://www.ebi.ac.uk/pride/archive/projects/PXD056286
- SQLite file format https://www.sqlite.org/fileformat.html

The manuscript remains suitable for assessment as exploratory public-data
research. This follow up strengthens traceability and interpretation. It does
not substitute for independent nerve validation or functional experiments.
