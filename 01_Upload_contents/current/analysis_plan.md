# Exploratory extension of the GBS and CIDP analysis

Plan recorded on 21 September 2026 before calculating the extension results.
The original manuscript results had already been examined. This is a post hoc
exploratory extension, not a prospectively registered confirmatory analysis.

## CSF evidence

Extract every CSF coefficient adjusted for age, sex and sample handling that is
printed in Tables 3, 4B, 5 and 6 of Kmezic et al. 2023
(doi:10.3389/fimmu.2023.1241199). Retain unavailable entries as missing.
Table 3 supplies IL8, which is not repeated in Table 4B. Include all proteins
with coefficients and standard errors in both independent GBS versus healthy
control cohorts in the two-cohort synthesis. Use the original fixed effect
inverse variance model, with BH correction across the entire expanded paired
protein family. Report DerSimonian and Laird random effects and modified
Hartung-Knapp intervals as sensitivities for every included protein, without
selecting a pooling model by significance. Preserve source-study P and FDR
values separately from normal-approximation calculations.

Present the source-study GBS versus CIDP and CIDP versus healthy control
coefficients separately, including all printed CSF proteins and IL8. These
contrasts share participants and are not independent validation cohorts or
additional inputs to the GBS meta-analysis. Retain the source-study adjusted
P values. Display normal-approximation intervals transparently, noting that
they need not reproduce the source small-sample regression tests.

## Nerve programs

Reconstruct pseudobulk expression from the deposited CellBender matrices and
original annotations for all 9 CIDP and 11 CIAP donors. Use whole-library
counts for CPM denominators. Require at least 20 nuclei per donor and cell
group, as in the original analysis. Standardize log2(CPM + 0.5) within each
cell group across eligible CIDP and CIAP donors with sample standard deviations.
All scores are equal-weight means of variable, measured genes. Missing genes
and zero-variance genes are reported and not replaced by zero.

Subsets are defined by biological function from the previously used panels,
before estimating the new contrasts. They do not purport to be validated
measures of permeability, transport flux, remyelination or functional repair.

- Endothelial junctions in ven_capEC2: CLDN5, OCLN, TJP1, CDH5.
- Endothelial transport-related transcripts in ven_capEC2: ABCB1, SLC1A1, MFSD2A.
- Endothelial inflammatory adhesion in ven_capEC2: ICAM1, VCAM1, SELE, SELP.
- Schwann myelin maintenance: MPZ, MBP, PRX, EGR2, PMP22, MAG.
- Schwann injury response: NGFR, ATF3, JUN, FOS, GDNF, RUNX2.

The two Schwann programs are examined in mySC, nmSC, and pooled repairSC plus
damageSC. SOX10 and S100B remain lineage markers rather than components of
either functional score. The structural and transport definitions are informed
by Palladino et al. 2017 (doi:10.1038/s41598-017-17475-y). The Schwann separation
is informed by Arthur-Farraj et al. 2012 (doi:10.1016/j.neuron.2012.06.021) and
Jessen and Mirsky 2016 (doi:10.1113/JP270874), already cited in the manuscript.

The nine new disease contrasts form one exploratory BH family. Also report
BH adjustment over these nine and all 80 archived nerve comparisons as a
conservative expanded-family sensitivity. Main P values use two-sided
Mann-Whitney tests with donors as independent units. Mean differences and
HC3 confidence intervals come from donor-level linear models, and their
distinct inferential methods are labelled. Age and centre adjusted HC3
models are prespecified here as sensitivities. For pooled repair and damage
Schwann cells, add the repairSC fraction to assess annotation composition.
Report all contrasts irrespective of direction or significance.

Relate the two new scores in pooled repair and damage Schwann cells to INCAT
within CIDP. Use Spearman correlations, BH over both new tests and, separately,
over those two plus the 15 original disability tests. These cross-sectional
correlations do not estimate recovery or treatment response.

## Reconstruction checks and attribution

Check original Fc gene-level mean differences and the archived Fc result
standardized within the 20 CIDP and CIAP donors against the raw reconstruction.
These are processing checks, not independent replication. Retain the original
80-comparison results and original mixed Schwann score as historical analyses.
Do not reinterpret changes in exploratory multiplicity as confirmation.

Normal BNB expression references and published tissue complement findings
provide anatomical or biological context. They do not independently replicate
the new Fc association. Any potential external disease dataset must be checked
for donor independence, cell type, comparator and endpoint before calling it
validation. No new wet-laboratory experiment is claimed.

## Subsequent small sample check

After examining the two new nine-donor clinical correlations, an exact
permutation sensitivity was added because the sample is small and INCAT has
tied values. It enumerates all 9! assignments and retains their multiplicities.
Both original asymptotic and exact new P values are retained. Main-text
interpretation uses the exact new tests. The expanded 17-test sensitivity
combines the 15 archived asymptotic P values with the two new exact P values
and is labelled accordingly. Leave-one-donor correlation ranges assess
influence within this dataset and are not independent replication.
