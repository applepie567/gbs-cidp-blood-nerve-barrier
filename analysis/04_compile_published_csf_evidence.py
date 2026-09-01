#!/usr/bin/env python3
"""Export the published GBS CSF evidence used in manuscript version 2.0.0.

This script preserves the inference level supplied by each publication. It does
not create participant-level observations from group-level published results.
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "source_data/Additional_file_1_source_data_v2.0.0.xlsx"
OUT = ROOT / "results/tables"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pxd = pd.read_excel(WORKBOOK, sheet_name="CSF_PXD002911")
    evidence = pd.read_excel(WORKBOOK, sheet_name="CSF_published_evidence")

    required = {
        "study", "resource", "cohort", "platform", "feature", "result",
        "program", "inference_level", "source",
    }
    missing = required.difference(evidence.columns)
    if missing:
        raise ValueError(f"Missing CSF evidence columns: {sorted(missing)}")
    if evidence["source"].isna().any() or evidence["inference_level"].isna().any():
        raise ValueError("Every CSF evidence row requires a source and inference level")

    pxd.to_csv(OUT / "csf_pxd002911_reported_results.csv", index=False)
    evidence.to_csv(OUT / "csf_published_evidence.csv", index=False)
    program_summary = (
        evidence.groupby("program", dropna=False)
        .agg(
            studies=("study", "nunique"),
            resources=("resource", lambda x: " | ".join(sorted(set(map(str, x))))),
            evidence_rows=("feature", "size"),
            inference_levels=("inference_level", lambda x: " | ".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )
    program_summary.to_csv(OUT / "csf_program_summary.csv", index=False)
    print(f"Exported {len(evidence)} published CSF evidence rows")


if __name__ == "__main__":
    main()
