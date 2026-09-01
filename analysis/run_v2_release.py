#!/usr/bin/env python3
"""Run the source-data compilation and release validation steps."""

from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "04_compile_published_csf_evidence.py",
    "05_compile_genetic_cell_localization.py",
    "06_build_cross_compartment_tables.py",
    "07_export_workbook_sheets.py",
]


def main() -> None:
    for step in STEPS:
        subprocess.run([sys.executable, str(ROOT / "analysis" / step)], check=True)
    subprocess.run([sys.executable, str(ROOT / "tests" / "validate_release.py")], check=True)


if __name__ == "__main__":
    main()
