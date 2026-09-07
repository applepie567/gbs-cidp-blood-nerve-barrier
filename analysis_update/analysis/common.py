#!/usr/bin/env python3
"""Shared prespecified panels and sample-level statistics for the upgrade.

The project intentionally uses the same compact pathway definitions across
platforms.  Tests are performed on biological samples (or matched patients),
never on individual cells or nuclei.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy import stats


GENE_MODULES = {
    "CXCL8_CXCR1_2": ["CXCL8", "CXCR1", "CXCR2"],
    "LIF_LIFR_IL6ST": [
        "LIF", "LIFR", "IL6", "IL6R", "IL6ST", "OSM", "OSMR",
        "STAT3", "SOCS3",
    ],
    "Complement": [
        "C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B", "C5",
        "CFB", "CFD", "CFH", "CFI", "SERPING1", "C3AR1", "C5AR1",
        "CR1", "CD55", "CD59",
    ],
    "Fc_receptor": [
        "FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B",
        "FCGRT", "FCER1G", "TYROBP", "SYK", "LYN", "HCK",
    ],
    "Transendothelial_migration": [
        "ITGAM", "ITGAL", "ITGA4", "ITGB1", "ITGB2", "SELL", "SELPLG",
        "CCR1", "CCR2", "CCR5", "CX3CR1", "ICAM1", "VCAM1", "SELE",
        "SELP", "PECAM1", "JAM2", "JAM3", "ALCAM", "MMP3", "MMP9",
        "CCL2", "CCL4", "CCL20",
    ],
    "Interferon_JAK_STAT": [
        "IFNB1", "IFNG", "STAT1", "STAT2", "IRF1", "IRF7", "ISG15",
        "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "JAK1", "JAK2",
    ],
    "Inflammatory_monocyte": [
        "CD14", "FCGR3A", "LYZ", "S100A8", "S100A9", "S100A12",
        "IL1B", "TNF", "NLRP3", "PYCARD", "OSM", "CD163", "C5AR1",
    ],
    "BNB_integrity": [
        "CLDN5", "OCLN", "TJP1", "CDH5", "PECAM1", "ABCB1", "MFSD2A",
        "SLC2A1", "GJA1", "PLVAP", "CAV1", "VWF", "EMCN", "ACKR1",
    ],
    "Schwann_myelin_repair": [
        "MPZ", "MBP", "PRX", "EGR2", "PMP22", "MAG", "NGFR", "ATF3",
        "JUN", "FOS", "GDNF", "RUNX2", "SOX10", "S100B",
    ],
}


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    indices = np.where(ok)[0]
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[indices[order]] = np.clip(adjusted, 0, 1)
    return out


def hedges_g(case: Iterable[float], reference: Iterable[float]) -> float:
    x = np.asarray(list(case), dtype=float)
    y = np.asarray(list(reference), dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float("nan")
    pooled = ((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2)
    if pooled <= 0 or not np.isfinite(pooled):
        return 0.0 if np.isclose(np.mean(x), np.mean(y)) else float("nan")
    correction = 1 - 3 / (4 * (nx + ny) - 9)
    return float(correction * (np.mean(x) - np.mean(y)) / math.sqrt(pooled))


def exact_permutation_p(case: Iterable[float], reference: Iterable[float]) -> float:
    """Two-sided exact label-permutation P for a mean difference.

    Enumeration is used only for the small cohorts in this project.  If the
    number of partitions exceeds 250,000, a deterministic 100,000 draw Monte
    Carlo approximation is used with the Phipson-Smyth +1 correction.
    """
    x = np.asarray(list(case), dtype=float)
    y = np.asarray(list(reference), dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    values = np.r_[x, y]
    n_case = len(x)
    if n_case == 0 or len(y) == 0:
        return float("nan")
    observed = abs(float(np.mean(x) - np.mean(y)))
    total = math.comb(len(values), n_case)
    tolerance = 1e-12
    if total <= 250_000:
        extreme = 0
        for idx in itertools.combinations(range(len(values)), n_case):
            mask = np.zeros(len(values), dtype=bool)
            mask[list(idx)] = True
            statistic = abs(float(np.mean(values[mask]) - np.mean(values[~mask])))
            extreme += statistic >= observed - tolerance
        return extreme / total
    rng = np.random.default_rng(20260821)
    extreme = 0
    for _ in range(100_000):
        perm = rng.permutation(len(values))
        statistic = abs(
            float(np.mean(values[perm[:n_case]]) - np.mean(values[perm[n_case:]]))
        )
        extreme += statistic >= observed - tolerance
    return (extreme + 1) / 100_001


def unpaired_effect(case: Iterable[float], reference: Iterable[float]) -> dict:
    x = np.asarray(list(case), dtype=float)
    y = np.asarray(list(reference), dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return {
            "n_case": len(x), "n_reference": len(y), "mean_case": np.nan,
            "mean_reference": np.nan, "delta": np.nan, "hedges_g": np.nan,
            "permutation_p": np.nan, "welch_p": np.nan,
        }
    return {
        "n_case": len(x),
        "n_reference": len(y),
        "mean_case": float(np.mean(x)),
        "mean_reference": float(np.mean(y)),
        "delta": float(np.mean(x) - np.mean(y)),
        "hedges_g": hedges_g(x, y),
        "permutation_p": exact_permutation_p(x, y),
        "welch_p": float(stats.ttest_ind(x, y, equal_var=False).pvalue),
    }


def paired_effect(case: Iterable[float], reference: Iterable[float]) -> dict:
    x = np.asarray(list(case), dtype=float)
    y = np.asarray(list(reference), dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    diff = x - y
    if len(diff) < 2:
        return {
            "n_pairs": len(diff), "mean_case": np.nan, "mean_reference": np.nan,
            "delta": np.nan, "dz": np.nan, "paired_t_p": np.nan,
            "wilcoxon_p": np.nan,
        }
    sd = float(np.std(diff, ddof=1))
    try:
        wilcoxon_p = float(stats.wilcoxon(diff, alternative="two-sided").pvalue)
    except ValueError:
        wilcoxon_p = 1.0
    return {
        "n_pairs": len(diff),
        "mean_case": float(np.mean(x)),
        "mean_reference": float(np.mean(y)),
        "delta": float(np.mean(diff)),
        "dz": float(np.mean(diff) / sd) if sd > 0 else 0.0,
        "paired_t_p": float(stats.ttest_rel(x, y).pvalue),
        "wilcoxon_p": wilcoxon_p,
    }


def standardized_module_scores(
    expression: pd.DataFrame,
    modules: dict[str, list[str]] | None = None,
    minimum_genes: int = 2,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Return sample-by-module scores from a gene-by-sample matrix.

    Each feature is standardized across biological samples before averaging;
    this keeps platform-specific abundance scales from dominating a module.
    """
    panels = modules or GENE_MODULES
    matrix = expression.apply(pd.to_numeric, errors="coerce")
    sd = matrix.std(axis=1, ddof=1).replace(0, np.nan)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(sd, axis=0)
    scores: dict[str, pd.Series] = {}
    available: dict[str, list[str]] = {}
    for module, genes in panels.items():
        present = [gene for gene in genes if gene in z.index]
        available[module] = present
        if len(present) >= minimum_genes:
            scores[module] = z.loc[present].mean(axis=0, skipna=True)
    return pd.DataFrame(scores), available

