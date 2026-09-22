#!/usr/bin/env python3
"""Targeted patient-level pseudobulk audit for GSE285983.

This is a feasibility analysis, not a definitive differential-expression model.
It aggregates raw counts by biological sample and curated nerve cell compartment,
then compares CIDP with CIAP and control nerves on log2(CPM + 0.5).
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[1]
H5_DIR = ROOT / "data/raw/GSE285983/h5"
METADATA = ROOT / "data/raw/GSE285983/GSE285983_metadata_all.csv.gz"
OUTPUT = ROOT / "results/gse285983_targeted_pseudobulk.json"
SAMPLE_EXPRESSION_OUTPUT = ROOT / "results/tables/gse285983_sample_targeted_expression.csv"
SAMPLE_MODULE_OUTPUT = ROOT / "results/tables/gse285983_sample_module_scores.csv"
SAMPLE_ABUNDANCE_OUTPUT = ROOT / "results/tables/gse285983_sample_cell_fractions.csv"


CELL_GROUPS = {
    "Macrophage": {"Macro1", "Macro2"},
    "Granulocyte": {"Granulo"},
    "T_NK": {"T_NK"},
    "B_cell": {"B"},
    "BNB_EC": {"ven_capEC2"},
    "Other_EC": {"artEC", "venEC", "ven_capEC1", "LEC"},
    "Pericyte": {"PC1", "PC2"},
    "Perineurium": {"periC1", "periC2", "periC3"},
    "Endoneurial_stroma": {"endoC"},
    "Myelinating_SC": {"mySC"},
    "Nonmyelinating_SC": {"nmSC"},
    "Repair_damage_SC": {"repairSC", "damageSC"},
}


GENE_PANELS = {
    "BNB_identity_integrity": [
        "ABCB1", "SLC1A1", "MFSD2A", "GJA1", "CLDN5", "OCLN", "TJP1",
        "CDH5", "PECAM1", "PLVAP", "CAV1", "VWF", "EMCN", "ACKR1",
    ],
    "Leukocyte_transmigration": [
        "ICAM1", "VCAM1", "SELE", "SELP", "JAM2", "JAM3", "ALCAM",
        "PECAM1", "CCL2", "CCR2", "CCL20", "CCR6", "MMP3", "MMP9",
        "ITGAM", "ITGAL", "ITGA4", "ITGB1", "ITGB2", "SELL", "SELPLG",
        "CCR1", "CCR5",
    ],
    "CXCL8_CXCR1_2": ["CXCL8", "CXCR1", "CXCR2"],
    "LIF_LIFR": [
        "LIF", "LIFR", "OSM", "OSMR", "IL6", "IL6R", "IL6ST", "STAT3", "SOCS3"
    ],
    "Complement": [
        "C1QA", "C1QB", "C1QC", "C3", "C5", "CFB", "CFD", "CFH", "CFI",
        "SERPING1", "C3AR1", "C5AR1", "CR1", "CR2", "CD55", "CD59",
    ],
    "Fc_receptor": [
        "FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B",
        "FCGRT", "FCER1G", "TYROBP", "SYK", "LYN", "HCK",
    ],
    "Macrophage_state": [
        "CD14", "CD163", "LYZ", "MS4A7", "TREM2", "SPP1", "APOE", "LPL",
        "FABP5", "GPNMB", "MRC1", "FOLR2", "CX3CR1",
    ],
    "Schwann_myelin_repair": [
        "MPZ", "MBP", "PRX", "EGR2", "PMP22", "MAG", "NGFR", "ATF3",
        "JUN", "FOS", "GDNF", "RUNX2", "SOX10", "S100B",
    ],
}


def bh_adjust(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg adjustment with monotonicity enforcement."""
    p = np.asarray(p_values, dtype=float)
    result = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not ok.any():
        return result.tolist()
    indices = np.where(ok)[0]
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result[indices[order]] = adjusted
    return result.tolist()


def read_metadata():
    sample_group = {}
    sample_incat = {}
    sample_metadata = {}
    barcode_cluster = defaultdict(dict)
    with gzip.open(METADATA, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            sample = row["sample"]
            barcode = row["barcode"]
            raw_barcode = barcode[len(sample) + 1 :]
            sample_group[sample] = row["level2"]
            sample_incat[sample] = row["incat"]
            sample_metadata[sample] = {
                "disease": row["level2"],
                "sex": row["sex"],
                "age": row["age"],
                "center": row["center"],
                "disease_duration_months": row["disease_duration_in_months"],
                "incat": row["incat"],
            }
            barcode_cluster[sample][raw_barcode] = row["cluster"]
    total_nuclei = {sample: len(items) for sample, items in barcode_cluster.items()}
    return sample_group, sample_incat, sample_metadata, barcode_cluster, total_nuclei


def aggregate_targeted_counts(barcode_cluster):
    target_genes = sorted({gene for genes in GENE_PANELS.values() for gene in genes})
    target_set = set(target_genes)
    aggregates = {}
    cell_counts = {}
    library_sizes = {}
    missing_genes = set(target_genes)

    for path in sorted(H5_DIR.glob("*.h5")):
        match = re.search(r"_(S\d+)_", path.name)
        if not match:
            raise ValueError(f"Cannot parse sample from {path.name}")
        sample = match.group(1)
        annotations = barcode_cluster[sample]

        with h5py.File(path, "r") as handle:
            matrix = handle["matrix"]
            gene_names = [x.decode() for x in matrix["features/name"][:]]
            gene_to_row = {gene: i for i, gene in enumerate(gene_names) if gene in target_set}
            missing_genes.difference_update(gene_to_row)
            selected_genes = sorted(gene_to_row)
            selected_rows = np.array([gene_to_row[g] for g in selected_genes], dtype=int)
            barcodes = [x.decode() for x in matrix["barcodes"][:]]

            data = matrix["data"][:]
            indices = matrix["indices"][:]
            indptr = matrix["indptr"][:]
            shape = tuple(int(x) for x in matrix["shape"][:])
            counts = sparse.csc_matrix((data, indices, indptr), shape=shape)
            targeted = counts[selected_rows, :]

            for group, clusters in CELL_GROUPS.items():
                columns = [
                    i
                    for i, barcode in enumerate(barcodes)
                    if annotations.get(barcode) in clusters
                ]
                cell_counts[f"{sample}|{group}"] = len(columns)
                if columns:
                    total = np.asarray(targeted[:, columns].sum(axis=1)).ravel()
                    library_sizes[f"{sample}|{group}"] = float(
                        counts[:, columns].sum()
                    )
                else:
                    total = np.zeros(len(selected_genes), dtype=float)
                    library_sizes[f"{sample}|{group}"] = 0.0
                aggregates[f"{sample}|{group}"] = {
                    gene: float(value) for gene, value in zip(selected_genes, total)
                }

    return target_genes, sorted(missing_genes), aggregates, cell_counts, library_sizes


def log_cpm(aggregate, genes, library_size):
    values = np.array([aggregate.get(g, 0.0) for g in genes], dtype=float)
    if library_size <= 0:
        return np.full(values.shape, np.nan)
    return np.log2(values / library_size * 1_000_000 + 0.5)


def targeted_comparisons(
    sample_group, genes, aggregates, cell_counts, library_sizes
):
    comparisons = [("CIDP", "CIAP"), ("CIDP", "CTRL")]
    rows = []

    for cell_group in CELL_GROUPS:
        profiles = {}
        for sample, disease in sample_group.items():
            key = f"{sample}|{cell_group}"
            if cell_counts.get(key, 0) < 20:
                continue
            profiles[sample] = log_cpm(
                aggregates[key], genes, library_sizes[key]
            )

        for case, reference in comparisons:
            case_samples = [s for s in profiles if sample_group[s] == case]
            ref_samples = [s for s in profiles if sample_group[s] == reference]
            if len(case_samples) < 2 or len(ref_samples) < 2:
                continue
            case_matrix = np.vstack([profiles[s] for s in case_samples])
            ref_matrix = np.vstack([profiles[s] for s in ref_samples])
            start = len(rows)
            for i, gene in enumerate(genes):
                x = case_matrix[:, i]
                y = ref_matrix[:, i]
                if np.allclose(np.r_[x, y], np.r_[x, y][0], equal_nan=False):
                    p_value = 1.0
                else:
                    p_value = float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
                rows.append(
                    {
                        "comparison": f"{case}_vs_{reference}",
                        "cell_group": cell_group,
                        "gene": gene,
                        "n_case": len(case_samples),
                        "n_reference": len(ref_samples),
                        "mean_log2cpm_case": float(np.mean(x)),
                        "mean_log2cpm_reference": float(np.mean(y)),
                        "delta_log2cpm": float(np.mean(x) - np.mean(y)),
                        "p_value": p_value,
                    }
                )
            adjusted = bh_adjust([row["p_value"] for row in rows[start:]])
            for row, q_value in zip(rows[start:], adjusted):
                row["fdr_within_celltype_comparison"] = q_value
    return rows


def expression_atlas(
    sample_group, genes, aggregates, cell_counts, library_sizes
):
    rows = []
    for disease in ["CIDP", "CIAP", "CTRL"]:
        for cell_group in CELL_GROUPS:
            sample_profiles = []
            contributing_samples = []
            for sample, group in sample_group.items():
                key = f"{sample}|{cell_group}"
                if group == disease and cell_counts.get(key, 0) >= 20:
                    sample_profiles.append(
                        log_cpm(aggregates[key], genes, library_sizes[key])
                    )
                    contributing_samples.append(sample)
            if not sample_profiles:
                continue
            matrix = np.vstack(sample_profiles)
            for index, gene in enumerate(genes):
                rows.append(
                    {
                        "disease": disease,
                        "cell_group": cell_group,
                        "gene": gene,
                        "n_samples": len(contributing_samples),
                        "mean_log2cpm": float(np.mean(matrix[:, index])),
                        "median_log2cpm": float(np.median(matrix[:, index])),
                    }
                )
    return rows


def pooled_expression_atlas(sample_group, genes, aggregates, cell_counts, library_sizes):
    """Descriptive pooled atlas for rare cell types; not used for inference."""
    rows = []
    for disease in ["CIDP", "CIAP", "CTRL"]:
        for cell_group in CELL_GROUPS:
            pooled_counts = {gene: 0.0 for gene in genes}
            pooled_library = 0.0
            pooled_nuclei = 0
            contributing_samples = 0
            for sample, group in sample_group.items():
                if group != disease:
                    continue
                key = f"{sample}|{cell_group}"
                if cell_counts.get(key, 0) <= 0:
                    continue
                contributing_samples += 1
                pooled_nuclei += cell_counts[key]
                pooled_library += library_sizes[key]
                for gene in genes:
                    pooled_counts[gene] += aggregates[key].get(gene, 0.0)
            if pooled_library <= 0:
                continue
            values = log_cpm(pooled_counts, genes, pooled_library)
            for gene, value in zip(genes, values):
                rows.append({
                    "disease": disease,
                    "cell_group": cell_group,
                    "gene": gene,
                    "pooled_log2cpm": float(value),
                    "pooled_nuclei": pooled_nuclei,
                    "contributing_samples": contributing_samples,
                })
    return rows


def module_comparisons(
    sample_group, genes, aggregates, cell_counts, library_sizes
):
    gene_index = {gene: i for i, gene in enumerate(genes)}
    comparisons = [("CIDP", "CIAP"), ("CIDP", "CTRL")]
    rows = []

    for cell_group in CELL_GROUPS:
        samples = []
        profiles = []
        for sample in sample_group:
            key = f"{sample}|{cell_group}"
            if cell_counts.get(key, 0) >= 20:
                samples.append(sample)
                profiles.append(
                    log_cpm(aggregates[key], genes, library_sizes[key])
                )
        if not profiles:
            continue
        matrix = np.vstack(profiles)
        means = np.mean(matrix, axis=0)
        stds = np.std(matrix, axis=0, ddof=1)
        stds[stds == 0] = np.nan
        z = (matrix - means) / stds

        for case, reference in comparisons:
            case_indices = [i for i, s in enumerate(samples) if sample_group[s] == case]
            ref_indices = [i for i, s in enumerate(samples) if sample_group[s] == reference]
            if len(case_indices) < 2 or len(ref_indices) < 2:
                continue
            start = len(rows)
            for panel, panel_genes in GENE_PANELS.items():
                indices = [gene_index[g] for g in panel_genes if g in gene_index]
                scores = np.nanmean(z[:, indices], axis=1)
                x = scores[case_indices]
                y = scores[ref_indices]
                p_value = float(
                    stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
                )
                rows.append(
                    {
                        "comparison": f"{case}_vs_{reference}",
                        "cell_group": cell_group,
                        "panel": panel,
                        "n_case": len(case_indices),
                        "n_reference": len(ref_indices),
                        "mean_z_case": float(np.mean(x)),
                        "mean_z_reference": float(np.mean(y)),
                        "delta_mean_z": float(np.mean(x) - np.mean(y)),
                        "p_value": p_value,
                    }
                )
            adjusted = bh_adjust([row["p_value"] for row in rows[start:]])
            for row, q_value in zip(rows[start:], adjusted):
                row["fdr_within_celltype_comparison"] = q_value
    return rows


def abundance_comparisons(sample_group, cell_counts, total_nuclei):
    comparisons = [("CIDP", "CIAP"), ("CIDP", "CTRL")]
    rows = []
    for case, reference in comparisons:
        start = len(rows)
        for cell_group in CELL_GROUPS:
            case_samples = [s for s, group in sample_group.items() if group == case]
            ref_samples = [s for s, group in sample_group.items() if group == reference]
            x = np.array(
                [cell_counts[f"{s}|{cell_group}"] / total_nuclei[s] for s in case_samples]
            )
            y = np.array(
                [cell_counts[f"{s}|{cell_group}"] / total_nuclei[s] for s in ref_samples]
            )
            p_value = float(
                stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
            )
            rows.append(
                {
                    "comparison": f"{case}_vs_{reference}",
                    "cell_group": cell_group,
                    "n_case": len(case_samples),
                    "n_reference": len(ref_samples),
                    "mean_fraction_case": float(np.mean(x)),
                    "mean_fraction_reference": float(np.mean(y)),
                    "difference_fraction": float(np.mean(x) - np.mean(y)),
                    "p_value": p_value,
                }
            )
        adjusted = bh_adjust([row["p_value"] for row in rows[start:]])
        for row, q_value in zip(rows[start:], adjusted):
            row["fdr_across_cell_groups"] = q_value
    return rows


def write_sample_level_tables(
    sample_group,
    sample_metadata,
    genes,
    aggregates,
    cell_counts,
    library_sizes,
    total_nuclei,
):
    """Write patient-level values used to construct manuscript figures."""
    expression_rows = []
    module_rows = []
    abundance_rows = []
    gene_index = {gene: i for i, gene in enumerate(genes)}

    for sample in sorted(sample_group):
        for cell_group in CELL_GROUPS:
            key = f"{sample}|{cell_group}"
            abundance_rows.append({
                "sample": sample,
                **sample_metadata[sample],
                "cell_group": cell_group,
                "n_nuclei": cell_counts.get(key, 0),
                "total_nuclei": total_nuclei[sample],
                "fraction": cell_counts.get(key, 0) / total_nuclei[sample],
            })
            if cell_counts.get(key, 0) < 20:
                continue
            values = log_cpm(aggregates[key], genes, library_sizes[key])
            for gene, value in zip(genes, values):
                expression_rows.append({
                    "sample": sample,
                    **sample_metadata[sample],
                    "cell_group": cell_group,
                    "gene": gene,
                    "log2_cpm_plus_0_5": float(value),
                    "n_nuclei": cell_counts[key],
                })

    # Module scores use z-standardized genes across all diseases within a cell group.
    for cell_group in CELL_GROUPS:
        samples = []
        profiles = []
        for sample in sorted(sample_group):
            key = f"{sample}|{cell_group}"
            if cell_counts.get(key, 0) >= 20:
                samples.append(sample)
                profiles.append(log_cpm(aggregates[key], genes, library_sizes[key]))
        if len(profiles) < 2:
            continue
        matrix = np.vstack(profiles)
        standard_deviation = np.std(matrix, axis=0, ddof=1)
        standard_deviation[standard_deviation == 0] = np.nan
        z = (matrix - np.mean(matrix, axis=0)) / standard_deviation
        for panel, panel_genes in GENE_PANELS.items():
            indices = [gene_index[gene] for gene in panel_genes if gene in gene_index]
            scores = np.nanmean(z[:, indices], axis=1)
            for sample, score in zip(samples, scores):
                module_rows.append({
                    "sample": sample,
                    **sample_metadata[sample],
                    "cell_group": cell_group,
                    "panel": panel,
                    "score_z": float(score),
                    "n_genes": len(indices),
                    "n_nuclei": cell_counts[f"{sample}|{cell_group}"],
                })

    SAMPLE_EXPRESSION_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    for path, rows in [
        (SAMPLE_EXPRESSION_OUTPUT, expression_rows),
        (SAMPLE_MODULE_OUTPUT, module_rows),
        (SAMPLE_ABUNDANCE_OUTPUT, abundance_rows),
    ]:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    return {
        "sample_expression_csv": str(SAMPLE_EXPRESSION_OUTPUT.relative_to(ROOT)),
        "sample_module_csv": str(SAMPLE_MODULE_OUTPUT.relative_to(ROOT)),
        "sample_abundance_csv": str(SAMPLE_ABUNDANCE_OUTPUT.relative_to(ROOT)),
    }


def main():
    (
        sample_group,
        sample_incat,
        sample_metadata,
        barcode_cluster,
        total_nuclei,
    ) = read_metadata()
    (
        genes,
        missing_genes,
        aggregates,
        cell_counts,
        library_sizes,
    ) = aggregate_targeted_counts(barcode_cluster)
    genes = [gene for gene in genes if gene not in missing_genes]
    comparisons = targeted_comparisons(
        sample_group, genes, aggregates, cell_counts, library_sizes
    )
    atlas = expression_atlas(
        sample_group, genes, aggregates, cell_counts, library_sizes
    )
    pooled_atlas = pooled_expression_atlas(
        sample_group, genes, aggregates, cell_counts, library_sizes
    )
    module_results = module_comparisons(
        sample_group, genes, aggregates, cell_counts, library_sizes
    )
    abundance_results = abundance_comparisons(
        sample_group, cell_counts, total_nuclei
    )
    sample_level_outputs = write_sample_level_tables(
        sample_group,
        sample_metadata,
        genes,
        aggregates,
        cell_counts,
        library_sizes,
        total_nuclei,
    )

    gene_to_panel = defaultdict(list)
    for panel, panel_genes in GENE_PANELS.items():
        for gene in panel_genes:
            gene_to_panel[gene].append(panel)
    for row in comparisons:
        row["panels"] = gene_to_panel[row["gene"]]
    for row in atlas:
        row["panels"] = gene_to_panel[row["gene"]]

    payload = {
        "dataset": "GSE285983",
        "analysis_level": "patient-level targeted pseudobulk",
        "normalization": "log2(CPM using all-gene library size + 0.5)",
        "minimum_nuclei_per_sample_cell_group": 20,
        "tests": "two-sided Mann-Whitney U; BH within cell group and comparison",
        "cautions": [
            "Exploratory feasibility analysis restricted to prespecified target genes.",
            "CIDP vs CTRL is confounded by sex and center in this cohort.",
            "CIDP vs CIAP is the preferred disease-control contrast but remains cross-sectional.",
            "Nuclear RNA can under-detect transient cytokine and chemokine transcripts.",
        ],
        "sample_groups": sample_group,
        "sample_incat": sample_incat,
        "sample_metadata": sample_metadata,
        "cell_groups": {k: sorted(v) for k, v in CELL_GROUPS.items()},
        "gene_panels": GENE_PANELS,
        "missing_genes": missing_genes,
        "cell_counts": cell_counts,
        "total_nuclei": total_nuclei,
        "library_sizes": library_sizes,
        "comparisons": comparisons,
        "module_comparisons": module_results,
        "abundance_comparisons": abundance_results,
        "expression_atlas": atlas,
        "pooled_expression_atlas": pooled_atlas,
        "sample_level_outputs": sample_level_outputs,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    print(f"Samples: {len(sample_group)}; comparisons: {len(comparisons)}")
    print(f"Missing genes: {', '.join(missing_genes) if missing_genes else 'none'}")


if __name__ == "__main__":
    main()
