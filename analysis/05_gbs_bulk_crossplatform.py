#!/usr/bin/env python3
"""Cross-platform bulk-transcriptome validation in human GBS blood.

GSE211225 contributes whole-blood RNA-seq from acute GBS, post-acute GBS,
and healthy controls. GSE31014 contributes an independent peripheral-
leukocyte microarray comparison. All tests use biological samples; module
definitions are identical to the single-cell analyses.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from common import GENE_MODULES, bh_adjust, standardized_module_scores, unpaired_effect  # noqa: E402


TABLE_ROOT = ROOT / "results/tables"
OUT_JSON = ROOT / "results/gbs_bulk_crossplatform.json"
ANALYSIS_MODULES = {
    key: value
    for key, value in GENE_MODULES.items()
    if key not in {"BNB_integrity", "Schwann_myelin_repair"}
}


def ensembl_symbol_map() -> dict[str, str]:
    path = ROOT / "data/raw/PRJNA1293757/HPBMC1/filtered_feature_bc_matrix.h5"
    with h5py.File(path, "r") as handle:
        ids = [value.decode().split(".")[0] for value in handle["matrix/features/id"][:]]
        names = [value.decode() for value in handle["matrix/features/name"][:]]
    return dict(zip(ids, names))


def median_ratio_normalize(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, str]:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    complete = (numeric > 0).all(axis=1)
    method = "DESeq-style median of ratios using genes nonzero in every sample"
    if int(complete.sum()) >= 100:
        geo = np.exp(np.log(numeric.loc[complete]).mean(axis=1))
        factors = numeric.loc[complete].div(geo, axis=0).median(axis=0)
    else:
        factors = numeric.sum(axis=0) / np.median(numeric.sum(axis=0))
        method = "total-count scaling fallback"
    factors = factors / np.exp(np.log(factors).mean())
    return numeric.div(factors, axis=1), factors, method


def parse_geo_metadata(path: Path) -> pd.DataFrame:
    wanted: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
                if row[0] in {"!Sample_geo_accession", "!Sample_title", "!Sample_description"}:
                    wanted[row[0]] = row[1:]
    out = pd.DataFrame({
        "gsm": wanted["!Sample_geo_accession"],
        "title": wanted["!Sample_title"],
        "description": wanted["!Sample_description"],
    })
    out["sample_number"] = out["description"].str.extract(r"/(\d+)\.bam")[0]
    out["condition"] = np.select(
        [out["title"].str.contains("Acute phase"), out["title"].str.contains("Post-acute")],
        ["Acute_GBS", "Postacute_GBS"],
        default="HC",
    )
    return out


def effect_tables(
    dataset: str,
    expression: pd.DataFrame,
    design: dict[str, str],
    comparisons: list[tuple[str, str, str]],
    modules: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    gene_rows: list[dict] = []
    module_rows: list[dict] = []
    score_rows: list[dict] = []
    scores, available = standardized_module_scores(expression, modules)
    for comparison, case_group, reference_group in comparisons:
        case = [sample for sample, group in design.items() if group == case_group]
        reference = [sample for sample, group in design.items() if group == reference_group]
        gene_start = len(gene_rows)
        for gene in expression.index:
            gene_rows.append({
                "dataset": dataset,
                "comparison": comparison,
                "gene": gene,
                **unpaired_effect(expression.loc[gene, case], expression.loc[gene, reference]),
                "modules": ";".join(name for name, genes in modules.items() if gene in genes),
            })
        for row, fdr in zip(gene_rows[gene_start:], bh_adjust(r["permutation_p"] for r in gene_rows[gene_start:])):
            row["targeted_fdr_within_comparison"] = fdr

        module_start = len(module_rows)
        for module in scores.columns:
            module_rows.append({
                "dataset": dataset,
                "comparison": comparison,
                "module": module,
                "n_genes": len(available[module]),
                "genes_available": ";".join(available[module]),
                **unpaired_effect(scores.loc[case, module], scores.loc[reference, module]),
            })
        for row, fdr in zip(module_rows[module_start:], bh_adjust(r["permutation_p"] for r in module_rows[module_start:])):
            row["module_fdr_within_comparison"] = fdr

    for sample in scores.index:
        for module in scores.columns:
            score_rows.append({
                "dataset": dataset,
                "sample": sample,
                "condition": design[sample],
                "module": module,
                "score_z": float(scores.loc[sample, module]),
            })
    return pd.DataFrame(gene_rows), pd.DataFrame(module_rows), pd.DataFrame(score_rows), available


def analyze_gse211225() -> dict:
    counts_path = ROOT / "data/raw/GSE211225/GSE211225_gene_counts_matrix_deseq.txt.gz"
    metadata = parse_geo_metadata(ROOT / "data/raw/GSE211225/GSE211225_series_matrix.txt.gz")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.columns = [re.search(r"/(\d+)\.bam", name).group(1) for name in counts.columns]
    mapping = ensembl_symbol_map()
    counts["gene"] = [mapping.get(str(index).split(".")[0], "") for index in counts.index]
    counts = counts.loc[counts["gene"].ne("")].groupby("gene", sort=False).sum(numeric_only=True)
    normalized, size_factors, normalization = median_ratio_normalize(counts)
    expression = np.log2(normalized + 0.5)
    target_genes = sorted({gene for genes in ANALYSIS_MODULES.values() for gene in genes})
    expression = expression.reindex([gene for gene in target_genes if gene in expression.index])
    sample_to_gsm = metadata.set_index("sample_number")["gsm"].to_dict()
    rename = {number: sample_to_gsm[number] for number in expression.columns}
    expression = expression.rename(columns=rename)
    design = metadata.set_index("gsm")["condition"].to_dict()
    comparisons = [
        ("Acute_GBS_vs_HC", "Acute_GBS", "HC"),
        ("Postacute_GBS_vs_HC", "Postacute_GBS", "HC"),
        ("Acute_GBS_vs_Postacute_GBS", "Acute_GBS", "Postacute_GBS"),
    ]
    genes, modules, scores, available = effect_tables(
        "GSE211225", expression, design, comparisons, ANALYSIS_MODULES
    )
    genes.to_csv(TABLE_ROOT / "gse211225_gene_effects.csv", index=False)
    modules.to_csv(TABLE_ROOT / "gse211225_module_effects.csv", index=False)
    scores.to_csv(TABLE_ROOT / "gse211225_sample_module_scores.csv", index=False)
    expression.rename_axis("gene").reset_index().to_csv(
        TABLE_ROOT / "gse211225_targeted_expression.csv", index=False
    )
    metadata.to_csv(TABLE_ROOT / "gse211225_sample_metadata.csv", index=False)
    return {
        "dataset": "GSE211225",
        "design_counts": metadata["condition"].value_counts().to_dict(),
        "normalization": normalization,
        "size_factors": {rename[k]: float(v) for k, v in size_factors.items()},
        "cross_sectional_not_paired": True,
        "module_availability": available,
        "gene_effects": genes.to_dict("records"),
        "module_effects": modules.to_dict("records"),
    }


def read_series_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))
    return pd.read_csv(io.StringIO("".join(lines[start:end])), sep="\t", index_col=0)


def read_gpl96_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("ID\tGene title\tGene symbol"))
    frame = pd.read_csv(io.StringIO("".join(lines[header:])), sep="\t", usecols=["ID", "Gene symbol"])
    frame["gene"] = frame["Gene symbol"].fillna("").str.split("///")
    frame = frame.explode("gene")
    frame["gene"] = frame["gene"].str.strip()
    return frame.loc[~frame["gene"].isin(["", "---"]), ["ID", "gene"]]


def analyze_gse31014() -> dict:
    probe = read_series_matrix(ROOT / "data/raw/GSE31014/GSE31014_series_matrix.txt.gz")
    annotation = read_gpl96_annotation(ROOT / "data/raw/GSE31014/GPL96.annot.gz")
    target_genes = sorted({gene for genes in ANALYSIS_MODULES.values() for gene in genes})
    annotation = annotation.loc[annotation["gene"].isin(target_genes)].copy()
    variance = probe.var(axis=1, ddof=1).rename("variance")
    annotation = annotation.join(variance, on="ID").sort_values(["gene", "variance"], ascending=[True, False])
    selected = annotation.drop_duplicates("gene").set_index("gene")["ID"]
    expression = probe.loc[selected.values].copy()
    expression.index = selected.index
    design = {sample: ("GBS" if index < 7 else "HC") for index, sample in enumerate(probe.columns)}
    genes, modules, scores, available = effect_tables(
        "GSE31014", expression, design, [("GBS_vs_HC", "GBS", "HC")], ANALYSIS_MODULES
    )
    selected.rename("selected_probe").reset_index().to_csv(TABLE_ROOT / "gse31014_probe_selection.csv", index=False)
    genes.to_csv(TABLE_ROOT / "gse31014_gene_effects.csv", index=False)
    modules.to_csv(TABLE_ROOT / "gse31014_module_effects.csv", index=False)
    scores.to_csv(TABLE_ROOT / "gse31014_sample_module_scores.csv", index=False)
    expression.rename_axis("gene").reset_index().to_csv(TABLE_ROOT / "gse31014_targeted_expression.csv", index=False)
    return {
        "dataset": "GSE31014",
        "design_counts": {"GBS": 7, "HC": 7},
        "input_scale": "GEO-deposited processed centered expression values",
        "probe_collapse": "highest intersample-variance probe per mapped gene",
        "module_availability": available,
        "gene_effects": genes.to_dict("records"),
        "module_effects": modules.to_dict("records"),
    }


def main() -> None:
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {"datasets": [analyze_gse211225(), analyze_gse31014()]}
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUT_JSON}")
    for dataset in payload["datasets"]:
        print(dataset["dataset"], dataset["design_counts"])


if __name__ == "__main__":
    main()
