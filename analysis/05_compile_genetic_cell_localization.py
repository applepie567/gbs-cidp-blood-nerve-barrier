#!/usr/bin/env python3
"""Export published CIDP genetic evidence and donor-resolved cell localization."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "source_data/Additional_file_1_source_data_v2.0.0.xlsx"
OUT = ROOT / "results/tables"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genetic = pd.read_excel(WORKBOOK, sheet_name="CIDP_genetic_evidence")
    donor = pd.read_excel(WORKBOOK, sheet_name="Genetic_donor_celltype")
    summary = pd.read_excel(WORKBOOK, sheet_name="Genetic_celltype_summary")
    contrast = pd.read_excel(WORKBOOK, sheet_name="Genetic_CIDP_vs_CIAP")
    tractability = pd.read_excel(WORKBOOK, sheet_name="OpenTargets_tractability")

    for name, frame in {
        "genetic": genetic,
        "donor": donor,
        "summary": summary,
        "contrast": contrast,
        "tractability": tractability,
    }.items():
        if "gene" not in frame.columns:
            raise ValueError(f"The {name} table does not contain a gene column")

    donor_counts = donor.groupby(["gene", "cell_group"], as_index=False).agg(
        donors=("sample", "nunique"),
        total_nuclei=("n_nuclei", "sum"),
        mean_expression=("mean_log2_cp10k_plus1", "mean"),
        mean_percent_expressing=("percent_expressing", "mean"),
    )
    donor_counts = donor_counts.sort_values(
        ["gene", "mean_percent_expressing", "mean_expression"],
        ascending=[True, False, False],
    )

    genetic.to_csv(OUT / "cidp_published_genetic_evidence.csv", index=False)
    donor.to_csv(OUT / "genetic_donor_celltype.csv", index=False)
    summary.to_csv(OUT / "genetic_celltype_summary.csv", index=False)
    contrast.to_csv(OUT / "genetic_cidp_vs_ciap.csv", index=False)
    tractability.to_csv(OUT / "opentargets_tractability.csv", index=False)
    donor_counts.to_csv(OUT / "genetic_celltype_priority_summary.csv", index=False)

    primary = set(genetic["gene"].dropna().astype(str))
    localized = set(donor["gene"].dropna().astype(str))
    absent = primary.difference(localized)
    if absent:
        raise ValueError(f"Published genetic genes missing from localization table: {sorted(absent)}")
    print(f"Localized {len(primary)} published CIDP genetic candidates")


if __name__ == "__main__":
    main()
