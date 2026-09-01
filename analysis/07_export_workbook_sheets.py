#!/usr/bin/env python3
"""Export every source-workbook sheet as a CSV file."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "source_data/Additional_file_1_source_data_v2.0.0.xlsx"
OUT = ROOT / "source_data/csv_by_sheet"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    book = pd.ExcelFile(WORKBOOK)
    for sheet in book.sheet_names:
        pd.read_excel(WORKBOOK, sheet_name=sheet).to_csv(OUT / f"{sheet}.csv", index=False)
    print(f"Exported {len(book.sheet_names)} workbook sheets")


if __name__ == "__main__":
    main()
