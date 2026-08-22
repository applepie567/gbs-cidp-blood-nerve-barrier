#!/usr/bin/env python3
"""Targeted patient-level reanalysis of early untreated GBS sorted-cell RNA-seq.

The GEO deposit GSE304871 provides normalized TPM-like expression matrices for
FACS-sorted CD11b+, CD4+, and CD8+ PBMC fractions.  Sample sizes are small, so
this script treats standardized mean differences and prespecified module scores
as the primary outputs; p values are descriptive and are not used to claim
genome-wide discovery.
"""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data/raw/GSE304871"
OUT_JSON = ROOT / "results/gse304871_targeted_sorted_bulk.json"
OUT_GENES = ROOT / "results/tables/gse304871_targeted_gene_effects.csv"
OUT_MODULES = ROOT / "results/tables/gse304871_module_effects.csv"
OUT_SCORES = ROOT / "results/tables/gse304871_sample_module_scores.csv"


SAMPLE_DESIGN = {
    "CD11b": {
        "356-1_S1": ("HC3", "HC"),
        "356-2_S2": ("HC2", "HC"),
        "356-3_S3": ("GBS3", "GBS"),
        "356-4_S4": ("HC1", "HC"),
        "356-5_S5": ("GBS1", "GBS"),
    },
    "CD4": {
        "356-6_S6": ("HC3", "HC"),
        "356-7_S7": ("HC2", "HC"),
        "356-8_S8": ("GBS3", "GBS"),
        "356-9_S9": ("HC1", "HC"),
        "356-10_S10": ("GBS1", "GBS"),
        "356-11_S11": ("GBS4", "GBS"),
    },
    "CD8": {
        "356-12_S12": ("HC3", "HC"),
        "356-13_S13": ("HC2", "HC"),
        "356-14_S14": ("GBS3", "GBS"),
        "356-15_S15": ("HC1", "HC"),
        "356-16_S16": ("GBS1", "GBS"),
    },
}


GENE_MODULES = {
    "CXCL8_axis": ["CXCL8", "CXCR1", "CXCR2"],
    "IL6_family_LIF_axis": ["LIF", "LIFR", "IL6", "IL6R", "IL6ST", "OSM", "OSMR", "STAT3", "SOCS3"],
    "Complement": ["C1QA", "C1QB", "C1QC", "C3", "C5", "CFB", "CFD", "CFH", "CFI", "SERPING1", "C3AR1", "C5AR1", "CR1", "CD55", "CD59"],
    "Fc_receptor": ["FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "FCGRT", "FCER1G", "TYROBP", "SYK", "LYN", "HCK"],
    "Transendothelial_migration": ["ITGAM", "ITGB2", "SELL", "SELPLG", "CCR2", "CCR1", "CCR5", "CX3CR1", "MMP9", "MMP3", "CCL2", "CCL4", "CCL20"],
    "Interferon_JAK_STAT": ["IFNB1", "IFNG", "STAT1", "STAT2", "IRF1", "IRF7", "ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "JAK1", "JAK2"],
    "Inflammatory_monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8", "S100A9", "S100A12", "IL1B", "TNF", "NLRP3", "OSM"],
}


def bh_adjust(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out.tolist()
    indices = np.where(ok)[0]
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out[indices[order]] = adjusted
    return out.tolist()


def hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    """Bias-corrected standardized mean difference (GBS minus healthy)."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float("nan")
    pooled_var = ((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2)
    if not np.isfinite(pooled_var) or pooled_var <= 0:
        return 0.0 if np.isclose(np.mean(x), np.mean(y)) else float("nan")
    d = (np.mean(x) - np.mean(y)) / math.sqrt(pooled_var)
    correction = 1 - 3 / (4 * (nx + ny) - 9)
    return float(correction * d)


def read_matrix(cell_type: str, gene_map: dict[str, str]) -> pd.DataFrame:
    path = DATA_DIR / f"GSE304871_{cell_type}_normalized_counts.csv.gz"
    with gzip.open(path, "rt") as handle:
        frame = pd.read_csv(handle, index_col=0)
    frame.index = frame.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    gene_names = frame.get("gene_name", pd.Series(index=frame.index, dtype=object)).copy()
    gene_names = gene_names.where(gene_names.notna(), frame.index.map(gene_map))
    sample_cols = list(SAMPLE_DESIGN[cell_type])
    expression = frame[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    expression.insert(0, "gene", gene_names.values)
    expression = expression.loc[expression["gene"].notna()]
    expression["gene"] = expression["gene"].astype(str)
    # In the unlikely event of duplicate symbols, retain their summed abundance.
    return expression.groupby("gene", sort=False)[sample_cols].sum()


def main() -> None:
    # CD4 has the complete symbol annotation and is used only as an identifier map.
    cd4_path = DATA_DIR / "GSE304871_CD4_normalized_counts.csv.gz"
    cd4 = pd.read_csv(
        cd4_path,
        compression="gzip",
        index_col=0,
        usecols=lambda column: column in {"Unnamed: 0", "gene_name"} or column == "",
    )
    cd4.index = cd4.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    gene_map = cd4["gene_name"].dropna().astype(str).to_dict()

    targeted = sorted({gene for genes in GENE_MODULES.values() for gene in genes})
    gene_rows: list[dict] = []
    module_rows: list[dict] = []
    score_rows: list[dict] = []
    missing: dict[str, list[str]] = {}

    for cell_type in ["CD11b", "CD4", "CD8"]:
        expression = read_matrix(cell_type, gene_map)
        present = [gene for gene in targeted if gene in expression.index]
        missing[cell_type] = sorted(set(targeted) - set(present))
        log_expr = np.log2(expression.loc[present] + 1.0)
        design = SAMPLE_DESIGN[cell_type]
        gbs_cols = [column for column, (_, group) in design.items() if group == "GBS"]
        hc_cols = [column for column, (_, group) in design.items() if group == "HC"]

        start = len(gene_rows)
        for gene in present:
            x = log_expr.loc[gene, gbs_cols].to_numpy(dtype=float)
            y = log_expr.loc[gene, hc_cols].to_numpy(dtype=float)
            p_value = float(stats.mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue)
            gene_rows.append({
                "cell_type": cell_type,
                "gene": gene,
                "n_gbs": len(x),
                "n_hc": len(y),
                "mean_log2_nTPM_plus1_gbs": float(np.mean(x)),
                "mean_log2_nTPM_plus1_hc": float(np.mean(y)),
                "delta_log2": float(np.mean(x) - np.mean(y)),
                "hedges_g": hedges_g(x, y),
                "p_value": p_value,
                "modules": ";".join(name for name, genes in GENE_MODULES.items() if gene in genes),
            })
        adjusted = bh_adjust([row["p_value"] for row in gene_rows[start:]])
        for row, q_value in zip(gene_rows[start:], adjusted):
            row["targeted_fdr_within_cell_type"] = q_value

        # Standardize genes across biological samples within each sorted fraction.
        standardized = log_expr.sub(log_expr.mean(axis=1), axis=0)
        standard_deviation = log_expr.std(axis=1, ddof=1).replace(0, np.nan)
        standardized = standardized.div(standard_deviation, axis=0)
        module_start = len(module_rows)
        for module, genes in GENE_MODULES.items():
            available = [gene for gene in genes if gene in standardized.index]
            if len(available) < 2:
                continue
            scores = standardized.loc[available].mean(axis=0, skipna=True)
            x = scores[gbs_cols].to_numpy(dtype=float)
            y = scores[hc_cols].to_numpy(dtype=float)
            p_value = float(stats.mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue)
            module_rows.append({
                "cell_type": cell_type,
                "module": module,
                "genes_available": ";".join(available),
                "n_genes": len(available),
                "n_gbs": len(x),
                "n_hc": len(y),
                "mean_z_gbs": float(np.mean(x)),
                "mean_z_hc": float(np.mean(y)),
                "delta_mean_z": float(np.mean(x) - np.mean(y)),
                "hedges_g": hedges_g(x, y),
                "p_value": p_value,
            })
            for column, score in scores.items():
                patient, group = design[column]
                score_rows.append({
                    "cell_type": cell_type,
                    "module": module,
                    "sample_column": column,
                    "patient": patient,
                    "group": group,
                    "score_z": float(score),
                })
        adjusted = bh_adjust([row["p_value"] for row in module_rows[module_start:]])
        for row, q_value in zip(module_rows[module_start:], adjusted):
            row["module_fdr_within_cell_type"] = q_value

    OUT_GENES.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(gene_rows).to_csv(OUT_GENES, index=False)
    pd.DataFrame(module_rows).to_csv(OUT_MODULES, index=False)
    pd.DataFrame(score_rows).to_csv(OUT_SCORES, index=False)

    payload = {
        "dataset": "GSE304871",
        "analysis_level": "patient-level targeted reanalysis of sorted-cell normalized expression",
        "transformation": "log2(nTPM + 1)",
        "primary_summary": "Hedges g and mean differences; exact Mann-Whitney tests are descriptive owing to very small n",
        "sample_design": SAMPLE_DESIGN,
        "gene_modules": GENE_MODULES,
        "missing_genes": missing,
        "gene_effects": gene_rows,
        "module_effects": module_rows,
        "cautions": [
            "Only two GBS donors are available for CD11b+ and CD8+ fractions and three for CD4+ cells.",
            "Normalized abundance rather than raw counts was deposited; no count-based dispersion model was fitted.",
            "The cohort represents early untreated AIDP-variant GBS and should not be generalized to all GBS subtypes.",
            "No recovery-phase sample is present in this GEO series.",
        ],
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUT_JSON}")
    print(f"Gene effects: {len(gene_rows)}; module effects: {len(module_rows)}")


if __name__ == "__main__":
    main()
