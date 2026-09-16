#!/usr/bin/env python3
"""Export published CIDP genetic evidence and donor-resolved cell localization."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
OUT = ROOT / "results/tables"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genetic = pd.read_csv(SRC / "CIDP_genetic_evidence.csv")
    donor = pd.read_csv(SRC / "Genetic_donor_celltype.csv")
    summary = pd.read_csv(SRC / "Genetic_celltype_summary.csv")
    contrast = pd.read_csv(SRC / "Genetic_CIDP_vs_CIAP.csv")

    for name, frame in {
        "genetic": genetic,
        "donor": donor,
        "summary": summary,
        "contrast": contrast,
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
    donor_counts.to_csv(OUT / "genetic_celltype_priority_summary.csv", index=False)

    primary = set(genetic["gene"].dropna().astype(str))
    localized = set(donor["gene"].dropna().astype(str))
    absent = primary.difference(localized)
    if absent:
        raise ValueError(f"Published genetic genes missing from localization table: {sorted(absent)}")
    print(f"Localized {len(primary)} published CIDP genetic candidates")


if __name__ == "__main__":
    main()
