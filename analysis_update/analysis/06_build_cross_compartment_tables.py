#!/usr/bin/env python3
"""Export the blood, CSF and peripheral-nerve evidence map."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
OUT = ROOT / "results/tables"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(SRC / "Cross_compartment_map.csv")
    blood = pd.read_csv(SRC / "Blood_cohort_effects.csv")

    required = {
        "Program", "Blood evidence", "CSF evidence",
        "Peripheral nerve evidence", "Interpretation",
    }
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(f"Missing cross-compartment columns: {sorted(missing)}")

    mapping.to_csv(OUT / "cross_compartment_evidence_map.csv", index=False)
    blood.to_csv(OUT / "blood_cohort_effects.csv", index=False)

    positive = (
        blood.assign(positive=blood["Positive direction"].astype(str).str.lower().eq("yes"))
        .groupby("Program", as_index=False)
        .agg(cohorts=("Cohort", "nunique"), positive_cohorts=("positive", "sum"))
    )
    positive.to_csv(OUT / "blood_direction_consistency.csv", index=False)
    print(f"Exported {len(mapping)} cross-compartment programs")


if __name__ == "__main__":
    main()
