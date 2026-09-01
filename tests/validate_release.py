#!/usr/bin/env python3
"""Validate the contents and analytical boundary of release v2.0.0."""

from pathlib import Path
import hashlib
import json
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "source_data/Additional_file_1_source_data_v2.0.0.xlsx"
FORBIDDEN = ["GSE304871", "26353948", "GBS-Proteomics", "longitudinal_proteomics"]
TEXT_EXTENSIONS = {".py", ".md", ".csv", ".json", ".yaml", ".yml", ".cff", ".txt"}


def scan_forbidden() -> list[dict]:
    hits = []
    scan_roots = [ROOT / "analysis", ROOT / "config", ROOT / "data", ROOT / "results", ROOT / "metadata"]
    for path in (p for base in scan_roots for p in base.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_EXTENSIONS:
            continue
        if path.name in {"FILES_TO_DELETE_FROM_V1.txt", "VALIDATION_REPORT.json"}:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for term in FORBIDDEN:
            if term.lower() in text.lower():
                hits.append({"file": str(path.relative_to(ROOT)), "term": term})
    return hits


def main() -> None:
    sheets = pd.ExcelFile(WORKBOOK).sheet_names
    required_sheets = {
        "Blood_cohort_effects", "CSF_PXD002911", "CSF_published_evidence",
        "CIDP_expression", "CIDP_module_effects", "Cross_compartment_map",
        "CIDP_genetic_evidence", "Genetic_donor_celltype",
        "Genetic_celltype_summary", "Genetic_CIDP_vs_CIAP",
    }
    missing = required_sheets.difference(sheets)
    figures = sorted((ROOT / "figures").glob("Figure_*.png"))
    tables = sorted((ROOT / "tables").glob("Table_*.csv"))
    dimensions = {}
    for figure in figures:
        with Image.open(figure) as image:
            dimensions[figure.name] = list(image.size)

    report = {
        "release": "2.0.0",
        "workbook_sheets": len(sheets),
        "missing_required_sheets": sorted(missing),
        "publication_figures": len(figures),
        "publication_tables": len(tables),
        "figure_dimensions": dimensions,
        "forbidden_content_hits": scan_forbidden(),
        "workbook_sha256": hashlib.sha256(WORKBOOK.read_bytes()).hexdigest(),
    }
    (ROOT / "metadata/VALIDATION_REPORT.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if missing or len(figures) != 5 or len(tables) != 3 or report["forbidden_content_hits"]:
        raise SystemExit(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
