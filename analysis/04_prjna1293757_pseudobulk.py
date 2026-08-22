#!/usr/bin/env python3
"""Independent acute GBS PBMC single-cell validation (PRJNA1293757).

Five public Cell Ranger matrices (three treatment-naive acute AIDP and two
matched healthy controls) are quality filtered and summarized as sample-level
pseudobulk profiles.  Cell-level P values are deliberately avoided.  A compact,
marker-based lineage classifier is used only to obtain a conservative monocyte
compartment for prespecified pathway scoring.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from common import GENE_MODULES, bh_adjust, standardized_module_scores, unpaired_effect  # noqa: E402


DATA_ROOT = ROOT / "data/raw/PRJNA1293757"
TABLE_ROOT = ROOT / "results/tables"
OUT_JSON = ROOT / "results/prjna1293757_pseudobulk.json"

SAMPLE_DESIGN = {
    "HPBMC1": "HC",
    "HPBMC2": "HC",
    "PBMC1": "Acute_GBS",
    "PBMC2": "Acute_GBS",
    "PBMC3": "Acute_GBS",
}

LINEAGE_MARKERS = {
    "Monocyte": ["LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LILRB1", "S100A8", "S100A9"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "LTB", "IL7R"],
    "B_cell": ["MS4A1", "CD79A", "CD37", "CD74", "HLA-DRA"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "Dendritic": ["FCER1A", "CD1C", "CLEC10A", "CST3"],
    "Platelet": ["PPBP", "PF4", "NRGN", "GNG11"],
}

ANALYSIS_MODULES = {
    key: value
    for key, value in GENE_MODULES.items()
    if key not in {"BNB_integrity", "Schwann_myelin_repair"}
}


def read_cellranger(path: Path):
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        names = np.array([item.decode() for item in matrix["features/name"][:]], dtype=object)
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(int(item) for item in matrix["shape"][:])
        counts = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return names, counts


def marker_annotation(
    names: np.ndarray,
    counts: sparse.csc_matrix,
    keep: np.ndarray,
    totals: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    name_to_row = {name: index for index, name in enumerate(names)}
    kept_columns = np.where(keep)[0]
    score_columns: dict[str, np.ndarray] = {}
    for lineage, markers in LINEAGE_MARKERS.items():
        rows = [name_to_row[gene] for gene in markers if gene in name_to_row]
        if not rows:
            score_columns[lineage] = np.full(len(kept_columns), -np.inf)
            continue
        selected = counts[rows, :][:, kept_columns].toarray().astype(float)
        selected = np.log1p(selected / totals[kept_columns] * 10_000)
        score_columns[lineage] = np.mean(selected, axis=0)
    score_frame = pd.DataFrame(score_columns)
    annotation = score_frame.idxmax(axis=1).to_numpy(dtype=object)
    return annotation, score_frame


def main() -> None:
    target_genes = sorted(
        {gene for genes in ANALYSIS_MODULES.values() for gene in genes}
        | {gene for genes in LINEAGE_MARKERS.values() for gene in genes}
    )
    sample_profiles: list[dict] = []
    qc_rows: list[dict] = []
    fraction_rows: list[dict] = []

    for sample, condition in SAMPLE_DESIGN.items():
        path = DATA_ROOT / sample / "filtered_feature_bc_matrix.h5"
        names, counts = read_cellranger(path)
        name_to_row = {name: index for index, name in enumerate(names)}
        totals = np.asarray(counts.sum(axis=0)).ravel().astype(float)
        detected = np.diff(counts.indptr)
        mito_rows = np.array([i for i, gene in enumerate(names) if str(gene).startswith("MT-")], dtype=int)
        mito = (
            np.asarray(counts[mito_rows, :].sum(axis=0)).ravel().astype(float)
            if len(mito_rows)
            else np.zeros(counts.shape[1], dtype=float)
        )
        mito_fraction = np.divide(
            mito, totals, out=np.zeros(mito.shape, dtype=float), where=totals > 0
        )
        keep = (detected >= 200) & (detected <= 6000) & (totals >= 500) & (mito_fraction < 0.20)
        annotation, _ = marker_annotation(names, counts, keep, totals)
        kept_columns = np.where(keep)[0]

        qc_rows.append({
            "dataset": "PRJNA1293757",
            "sample": sample,
            "condition": condition,
            "cellranger_filtered_cells": int(counts.shape[1]),
            "qc_retained_cells": int(keep.sum()),
            "median_umi_retained": float(np.median(totals[keep])),
            "median_genes_retained": float(np.median(detected[keep])),
            "median_mito_fraction_retained": float(np.median(mito_fraction[keep])),
        })

        targeted_present = [gene for gene in target_genes if gene in name_to_row]
        targeted_rows = np.array([name_to_row[gene] for gene in targeted_present], dtype=int)
        for compartment, column_mask in {
            "All_PBMC": np.ones(len(kept_columns), dtype=bool),
            "Marker_defined_monocyte": annotation == "Monocyte",
        }.items():
            selected_columns = kept_columns[column_mask]
            library_size = float(counts[:, selected_columns].sum())
            target_sum = np.asarray(counts[targeted_rows, :][:, selected_columns].sum(axis=1)).ravel()
            log_cpm = np.log2(target_sum / library_size * 1_000_000 + 0.5)
            for gene, value in zip(targeted_present, log_cpm):
                sample_profiles.append({
                    "dataset": "PRJNA1293757",
                    "sample": sample,
                    "condition": condition,
                    "compartment": compartment,
                    "gene": gene,
                    "log2_cpm_plus_0_5": float(value),
                    "n_cells": int(len(selected_columns)),
                    "library_size": library_size,
                })

        for lineage in LINEAGE_MARKERS:
            fraction_rows.append({
                "dataset": "PRJNA1293757",
                "sample": sample,
                "condition": condition,
                "lineage": lineage,
                "n_cells": int(np.sum(annotation == lineage)),
                "fraction_qc_cells": float(np.mean(annotation == lineage)),
            })

    profile_frame = pd.DataFrame(sample_profiles)
    qc_frame = pd.DataFrame(qc_rows)
    fraction_frame = pd.DataFrame(fraction_rows)
    gene_effect_rows: list[dict] = []
    module_effect_rows: list[dict] = []
    module_score_rows: list[dict] = []
    module_availability: dict[str, dict[str, list[str]]] = {}

    for compartment, subset in profile_frame.groupby("compartment"):
        expression = subset.pivot(index="gene", columns="sample", values="log2_cpm_plus_0_5")
        design = subset.drop_duplicates("sample").set_index("sample")["condition"].to_dict()
        acute = [sample for sample, group in design.items() if group == "Acute_GBS"]
        controls = [sample for sample, group in design.items() if group == "HC"]
        start = len(gene_effect_rows)
        for gene in expression.index:
            effect = unpaired_effect(expression.loc[gene, acute], expression.loc[gene, controls])
            gene_effect_rows.append({
                "dataset": "PRJNA1293757",
                "comparison": "Acute_GBS_vs_HC",
                "compartment": compartment,
                "gene": gene,
                **effect,
                "modules": ";".join(k for k, genes in ANALYSIS_MODULES.items() if gene in genes),
            })
        q = bh_adjust(row["permutation_p"] for row in gene_effect_rows[start:])
        for row, value in zip(gene_effect_rows[start:], q):
            row["targeted_fdr_within_compartment"] = value

        scores, available = standardized_module_scores(expression, ANALYSIS_MODULES)
        module_availability[compartment] = available
        start = len(module_effect_rows)
        for module in scores.columns:
            effect = unpaired_effect(scores.loc[acute, module], scores.loc[controls, module])
            module_effect_rows.append({
                "dataset": "PRJNA1293757",
                "comparison": "Acute_GBS_vs_HC",
                "compartment": compartment,
                "module": module,
                "n_genes": len(available[module]),
                "genes_available": ";".join(available[module]),
                **effect,
            })
            for sample, score in scores[module].items():
                module_score_rows.append({
                    "dataset": "PRJNA1293757",
                    "sample": sample,
                    "condition": design[sample],
                    "compartment": compartment,
                    "module": module,
                    "score_z": float(score),
                })
        q = bh_adjust(row["permutation_p"] for row in module_effect_rows[start:])
        for row, value in zip(module_effect_rows[start:], q):
            row["module_fdr_within_compartment"] = value

    fraction_effect_rows = []
    for lineage, subset in fraction_frame.groupby("lineage"):
        fraction_effect_rows.append({
            "dataset": "PRJNA1293757",
            "comparison": "Acute_GBS_vs_HC",
            "lineage": lineage,
            **unpaired_effect(
                subset.loc[subset["condition"] == "Acute_GBS", "fraction_qc_cells"],
                subset.loc[subset["condition"] == "HC", "fraction_qc_cells"],
            ),
        })

    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    profile_frame.to_csv(TABLE_ROOT / "prjna1293757_sample_targeted_expression.csv", index=False)
    qc_frame.to_csv(TABLE_ROOT / "prjna1293757_qc_summary.csv", index=False)
    fraction_frame.to_csv(TABLE_ROOT / "prjna1293757_sample_cell_fractions.csv", index=False)
    pd.DataFrame(fraction_effect_rows).to_csv(TABLE_ROOT / "prjna1293757_cell_fraction_effects.csv", index=False)
    pd.DataFrame(gene_effect_rows).to_csv(TABLE_ROOT / "prjna1293757_gene_effects.csv", index=False)
    pd.DataFrame(module_effect_rows).to_csv(TABLE_ROOT / "prjna1293757_module_effects.csv", index=False)
    pd.DataFrame(module_score_rows).to_csv(TABLE_ROOT / "prjna1293757_sample_module_scores.csv", index=False)

    payload = {
        "dataset": "PRJNA1293757",
        "design": SAMPLE_DESIGN,
        "analysis_level": "biological-sample pseudobulk; no cell-level hypothesis tests",
        "qc": qc_rows,
        "lineage_annotation": "highest mean log-normalized expression among prespecified canonical lineage marker panels",
        "module_availability": module_availability,
        "gene_effects": gene_effect_rows,
        "module_effects": module_effect_rows,
        "cell_fraction_effects": fraction_effect_rows,
        "cautions": [
            "Only three acute AIDP cases and two controls are available; exact label-permutation P values are coarse because only 10 case/control partitions exist.",
            "Marker-defined monocytes provide a conservative validation compartment and are not a substitute for the authors' full reference-based annotation.",
            "The cohort is independent of GSE304871/GSE304872 and represents treatment-naive acute AIDP.",
        ],
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUT_JSON}")
    print(qc_frame[["sample", "condition", "qc_retained_cells"]].to_string(index=False))


if __name__ == "__main__":
    main()
