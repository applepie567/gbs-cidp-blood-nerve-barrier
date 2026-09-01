#!/usr/bin/env python3
"""Export the blood, CSF and peripheral-nerve evidence map."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "source_data/Additional_file_1_source_data_v2.0.0.xlsx"
OUT = ROOT / "results/tables"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_excel(WORKBOOK, sheet_name="Cross_compartment_map")
    blood = pd.read_excel(WORKBOOK, sheet_name="Blood_cohort_effects")
    figure_index = pd.read_excel(WORKBOOK, sheet_name="Figure_table_index")

    required = {
        "Program", "Blood evidence", "CSF evidence",
        "Peripheral-nerve evidence", "Manuscript interpretation",
    }
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(f"Missing cross-compartment columns: {sorted(missing)}")

    mapping.to_csv(OUT / "cross_compartment_evidence_map.csv", index=False)
    blood.to_csv(OUT / "blood_cohort_effects.csv", index=False)
    figure_index.to_csv(OUT / "figure_table_source_index.csv", index=False)

    positive = (
        blood.assign(positive=blood["Positive direction"].astype(str).str.lower().eq("yes"))
        .groupby("Program", as_index=False)
        .agg(cohorts=("Cohort", "nunique"), positive_cohorts=("positive", "sum"))
    )
    positive.to_csv(OUT / "blood_direction_consistency.csv", index=False)
    print(f"Exported {len(mapping)} cross-compartment programs")


if __name__ == "__main__":
    main()
