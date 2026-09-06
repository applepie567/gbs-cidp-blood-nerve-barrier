#!/usr/bin/env python3
"""Verify public aggregate results, or recalculate with explicit local inputs."""

from pathlib import Path
import subprocess
import sys
import argparse

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "04_compile_published_csf_evidence.py",
    "05_compile_genetic_cell_localization.py",
    "06_build_cross_compartment_tables.py",
    "07_export_workbook_sheets.py",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-local-donor-tables', action='store_true',
                        help='Recalculate donor analyses after reconstructing local inputs.')
    args = parser.parse_args()
    steps = ['13_verify_public_summaries.py']
    if args.with_local_donor_tables:
        required = ['source_data/Genetic_donor_celltype.csv',
                    'source_data/CIDP_donor_heterogeneity.csv',
                    'results/tables/gse285983_sample_module_scores.csv']
        missing = [p for p in required if not (ROOT/p).is_file()]
        if missing:
            parser.error('Required local donor inputs are missing: ' + ', '.join(missing)
                         + '. See metadata/data_reconstruction.csv and Additional file 2.')
        steps = STEPS
    for step in steps:
        subprocess.run([sys.executable, str(ROOT / "analysis" / step)], check=True)
    subprocess.run([sys.executable, str(ROOT / "tests" / "validate_release.py")], check=True)


if __name__ == "__main__":
    main()
