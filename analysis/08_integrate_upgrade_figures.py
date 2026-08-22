#!/usr/bin/env python3
"""Integrate upgraded cohorts and generate manuscript figures/tables."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "results/tables"
FIG = ROOT / "submission/figures"
OUT_JSON = ROOT / "results/manuscript_stats_upgraded.json"

COLORS = {
    "ink": "#17212B", "blue": "#2474A6", "light_blue": "#DCECF5",
    "red": "#C64B47", "orange": "#E18A3B", "green": "#3B9672",
    "purple": "#795C9D", "gray": "#6B7280", "light": "#F2F5F7",
    "line": "#D8DEE3", "gold": "#B98B2F",
}

MODULE_LABELS = {
    "CXCL8_axis": "CXCL8-CXCR1/2",
    "CXCL8_CXCR1_2": "CXCL8-CXCR1/2",
    "IL6_family_LIF_axis": "LIF/OSM-gp130",
    "LIF_LIFR_IL6ST": "LIF/OSM-gp130",
    "LIF_LIFR": "LIF/OSM-gp130",
    "Complement": "Complement",
    "Fc_receptor": "Fc receptor",
    "Transendothelial_migration": "Transmigration",
    "Leukocyte_transmigration": "Transmigration",
    "Interferon_JAK_STAT": "IFN-JAK-STAT",
    "Inflammatory_monocyte": "Inflammatory monocyte",
    "Macrophage_state": "Macrophage state",
    "BNB_identity_integrity": "BNB integrity",
    "Schwann_myelin_repair": "Schwann/myelin repair",
}
CORE = [
    "CXCL8-CXCR1/2", "LIF/OSM-gp130", "Complement", "Fc receptor",
    "Transmigration", "IFN-JAK-STAT", "Inflammatory monocyte",
]


def style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.2,
        "axes.titlesize": 9.5, "axes.labelsize": 8.2,
        "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
        "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    sns.set_theme(style="whitegrid", rc={"grid.color": "#E6EAED", "grid.linewidth": 0.5})


def panel(ax, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", color=COLORS["ink"])


def save(fig: plt.Figure, stem: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def module_effect_matrix() -> pd.DataFrame:
    rows: list[dict] = []
    g304 = pd.read_csv(TABLE / "gse304871_module_effects.csv")
    for row in g304.loc[g304["cell_type"].eq("CD11b")].itertuples(index=False):
        rows.append({"cohort": "GSE304871 CD11b+", "contrast": "acute GBS vs HC",
                     "module": MODULE_LABELS[row.module], "effect": row.hedges_g,
                     "p": row.p_value, "fdr": row.module_fdr_within_cell_type,
                     "effect_metric": "Hedges g", "evidence": "direct human transcriptome"})

    prj = pd.read_csv(TABLE / "prjna1293757_module_effects.csv")
    prj = prj.loc[prj["compartment"].eq("Marker_defined_monocyte")]
    for row in prj.itertuples(index=False):
        rows.append({"cohort": "PRJNA1293757 monocyte", "contrast": "acute GBS vs HC",
                     "module": MODULE_LABELS[row.module], "effect": row.hedges_g,
                     "p": row.permutation_p, "fdr": row.module_fdr_within_compartment,
                     "effect_metric": "Hedges g", "evidence": "direct human scRNA pseudobulk"})

    g211 = pd.read_csv(TABLE / "gse211225_module_effects.csv")
    for comparison, label in {
        "Acute_GBS_vs_HC": "acute GBS vs HC",
        "Acute_GBS_vs_Postacute_GBS": "acute vs post-acute GBS",
    }.items():
        for row in g211.loc[g211["comparison"].eq(comparison)].itertuples(index=False):
            rows.append({"cohort": "GSE211225 whole blood", "contrast": label,
                         "module": MODULE_LABELS[row.module], "effect": row.hedges_g,
                         "p": row.permutation_p, "fdr": row.module_fdr_within_comparison,
                         "effect_metric": "Hedges g", "evidence": "direct human transcriptome"})

    g310 = pd.read_csv(TABLE / "gse31014_module_effects.csv")
    for row in g310.itertuples(index=False):
        rows.append({"cohort": "GSE31014 leukocytes", "contrast": "GBS vs HC",
                     "module": MODULE_LABELS[row.module], "effect": row.hedges_g,
                     "p": row.permutation_p, "fdr": row.module_fdr_within_comparison,
                     "effect_metric": "Hedges g", "evidence": "direct human microarray"})

    prot = pd.read_csv(TABLE / "gbs_proteomics_module_effects.csv")
    for comparison, label in {
        "Acute_GBS_vs_HC": "acute GBS vs HC",
        "Acute_GBS_vs_Recovery_1y": "acute vs paired 1-y recovery",
    }.items():
        for row in prot.loc[prot["comparison"].eq(comparison)].itertuples(index=False):
            paired = comparison.endswith("Recovery_1y")
            rows.append({"cohort": "GBS-Proteomics plasma", "contrast": label,
                         "module": MODULE_LABELS[row.module],
                         "effect": row.dz if paired else row.hedges_g,
                         "p": row.paired_t_p if paired else row.welch_p,
                         "fdr": row.module_fdr_within_comparison,
                         "effect_metric": "paired dz" if paired else "Hedges g",
                         "evidence": "direct human proteome; preprint"})
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLE / "upgraded_module_effect_matrix.csv", index=False)
    return frame


def key_stats(effect: pd.DataFrame) -> dict:
    g211 = pd.read_csv(TABLE / "gse211225_module_effects.csv")
    prot = pd.read_csv(TABLE / "gbs_proteomics_module_effects.csv")
    cidp = pd.read_csv(TABLE / "gse285983_all_module_effects.csv")
    bnb = pd.read_csv(TABLE / "gse107574_bnb_target_summary.csv")
    prj_qc = pd.read_csv(TABLE / "prjna1293757_qc_summary.csv")
    oep = pd.read_csv(TABLE / "oep002315_oep002701_selected_author_gsea.csv")

    acute = g211.loc[g211["comparison"].eq("Acute_GBS_vs_HC")].set_index("module")
    acute_post = g211.loc[g211["comparison"].eq("Acute_GBS_vs_Postacute_GBS")].set_index("module")
    paired = prot.loc[prot["comparison"].eq("Acute_GBS_vs_Recovery_1y")].set_index("module")
    fc = cidp.loc[
        cidp["comparison"].eq("CIDP_vs_CIAP")
        & cidp["cell_group"].eq("Macrophage")
        & cidp["panel"].eq("Fc_receptor")
    ].iloc[0]
    bnb_key = {}
    for gene in ["CLDN5", "OCLN", "TJP1", "CDH5", "VWF", "LIFR", "IL6ST", "ICAM1", "VCAM1", "C3"]:
        subset = bnb.loc[bnb["gene"].eq(gene)].set_index("preparation")
        bnb_key[gene] = {
            "LCM_median_fpkm": float(subset.loc["LCM_endoneurial_microvessel", "median_fpkm"]),
            "cultured_EC_median_fpkm": float(subset.loc["cultured_endoneurial_endothelial_cell", "median_fpkm"]),
        }
    payload = {
        "direct_human_resources": ["GSE304871", "PRJNA1293757", "GSE211225", "GSE31014", "GBS-Proteomics", "GSE285983", "GSE107574"],
        "direct_human_resource_count": 7,
        "gse211225_design": {"acute": 6, "postacute": 10, "healthy": 6, "paired": False},
        "gse211225_acute_vs_hc": {
            module: {"delta": float(acute.loc[module, "delta"]),
                     "hedges_g": float(acute.loc[module, "hedges_g"]),
                     "permutation_p": float(acute.loc[module, "permutation_p"]),
                     "fdr": float(acute.loc[module, "module_fdr_within_comparison"])}
            for module in acute.index
        },
        "gse211225_acute_vs_postacute": {
            module: {"delta": float(acute_post.loc[module, "delta"]),
                     "hedges_g": float(acute_post.loc[module, "hedges_g"]),
                     "permutation_p": float(acute_post.loc[module, "permutation_p"]),
                     "fdr": float(acute_post.loc[module, "module_fdr_within_comparison"])}
            for module in acute_post.index
        },
        "prjna1293757_qc_retained_cells": int(prj_qc["qc_retained_cells"].sum()),
        "gbs_proteomics": {
            "design": {"acute": 20, "paired_recovery_1y": 20, "healthy": 15},
            "paired_effects": {
                module: {"delta": float(paired.loc[module, "delta"]),
                         "dz": float(paired.loc[module, "dz"]),
                         "paired_t_p": float(paired.loc[module, "paired_t_p"]),
                         "fdr": float(paired.loc[module, "module_fdr_within_comparison"])}
                for module in paired.index
            },
        },
        "cidp_macrophage_fc_receptor": {
            "delta_z": float(fc["delta_mean_z"]), "p": float(fc["p_value"]),
            "fdr": float(fc["fdr_within_celltype_comparison"]),
        },
        "normal_bnb_anchor": bnb_key,
        "oep_cd16_complement_gsea": oep.loc[
            oep["cell_subset"].eq("CD16+ monocytes KEGG")
            & oep["description"].eq("Complement and coagulation cascades")
        ].iloc[0].to_dict(),
        "effect_matrix_rows": effect.to_dict("records"),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return payload


def figure1() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.set_xlim(0, 12); ax.set_ylim(0, 10); ax.axis("off")
    panel(ax, "a")
    ax.text(0, 10.15, "Upgraded compartment- and phase-aware evidence architecture",
            fontsize=10.5, fontweight="bold", color=COLORS["ink"])
    columns = [
        (0.2, 2.2, "Acute GBS\nblood", COLORS["red"]),
        (3.0, 2.25, "GBS post-acute /\n1-y recovery", COLORS["orange"]),
        (6.05, 2.25, "BNB interface", COLORS["blue"]),
        (9.0, 2.7, "CIDP peripheral\nnerve", COLORS["green"]),
    ]
    for x, width, title, color in columns:
        ax.add_patch(FancyBboxPatch((x, 2.0), width, 6.8, boxstyle="round,pad=0.08",
                                    facecolor="white", edgecolor=color, lw=1.5))
        ax.add_patch(Rectangle((x, 7.75), width, 1.05, facecolor=color, edgecolor=color))
        ax.text(x + width/2, 8.28, title, ha="center", va="center", color="white",
                fontsize=8.2, fontweight="bold")
    items = [
        (0.4, 7.15, "GSE211225 whole blood\n6 acute vs 6 HC", "direct"),
        (0.4, 6.05, "PRJNA1293757 scRNA\n3 acute vs 2 HC", "direct"),
        (0.4, 4.95, "GSE304871 CD11b+\n2 acute vs 3 HC", "direct"),
        (0.4, 3.85, "GSE31014 leukocytes\n7 GBS vs 7 HC", "direct"),
        (3.25, 7.15, "GSE211225 post-acute\n10 (cross-sectional)", "direct"),
        (3.25, 5.85, "GBS-Proteomics\n20 acute / recovery pairs", "direct preprint"),
        (6.3, 7.15, "GSE107574 normal BNB\n2 EC + 4 microvessels", "identity anchor"),
        (6.3, 5.75, "GSE285983 CIDP BNB EC\ndonor pseudobulk", "direct"),
        (9.25, 7.15, "GSE285983\n9 CIDP; 37 nerves", "direct"),
        (9.25, 5.75, "2026 stable-treated CIDP\nstudy; publication-level", "external"),
        (9.25, 4.35, "IVIg response cohorts\naccess / heterogeneity gate", "guardrail"),
    ]
    tag_color = {"direct": COLORS["green"], "direct preprint": COLORS["gold"],
                 "identity anchor": COLORS["blue"], "external": COLORS["orange"],
                 "guardrail": COLORS["gray"]}
    for x, y, text, tag in items:
        ax.text(x, y, text, va="top", fontsize=6.9, color=COLORS["ink"])
        ax.text(x, y-0.63, tag, va="top", fontsize=5.8, fontweight="bold", color=tag_color[tag])
    for start, end in [
        ((2.4, 1.7), (3.0, 1.7)), ((5.25, 1.7), (6.05, 1.7)),
        ((8.3, 1.7), (9.0, 1.7)),
    ]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11,
                                     color=COLORS["gray"], lw=1.2))
    ax.add_patch(FancyBboxPatch((0.35, 0.25), 11.3, 1.05, boxstyle="round,pad=0.05",
                                fc=COLORS["light"], ec=COLORS["line"], lw=1))
    ax.text(6, 0.77, "Anatomical guardrail: peripheral endoneurial microvascular endothelium = BNB, not BBB;\n"
            "the perineurium remains a separate PNS diffusion barrier.", ha="center", va="center",
            fontsize=6.8, fontweight="bold", color=COLORS["red"])
    save(fig, "Figure_1_study_design")


def figure2(effect: pd.DataFrame) -> None:
    row_order = [
        ("GSE304871 CD11b+", "acute GBS vs HC"),
        ("PRJNA1293757 monocyte", "acute GBS vs HC"),
        ("GSE211225 whole blood", "acute GBS vs HC"),
        ("GSE211225 whole blood", "acute vs post-acute GBS"),
        ("GSE31014 leukocytes", "GBS vs HC"),
        ("GBS-Proteomics plasma", "acute GBS vs HC"),
        ("GBS-Proteomics plasma", "acute vs paired 1-y recovery"),
    ]
    labels = [f"{a}\n{b}" for a, b in row_order]
    effect["row"] = effect["cohort"] + "\n" + effect["contrast"]
    matrix = effect.pivot_table(index="row", columns="module", values="effect", aggfunc="first").reindex(index=labels, columns=CORE)
    fdr = effect.pivot_table(index="row", columns="module", values="fdr", aggfunc="first").reindex(index=labels, columns=CORE)
    pvals = effect.pivot_table(index="row", columns="module", values="p", aggfunc="first").reindex(index=labels, columns=CORE)
    fig, ax = plt.subplots(figsize=(7.2, 6.15))
    fig.subplots_adjust(left=.29, right=.91, top=.91, bottom=.28)
    panel(ax, "a")
    sns.heatmap(matrix, cmap="vlag", center=0, vmin=-2.5, vmax=2.5, linewidths=.6,
                linecolor="white", mask=matrix.isna(), cbar_kws={"label": "Standardized effect", "shrink": .72}, ax=ax)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if not np.isfinite(value):
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#E5E7EB", edgecolor="white", lw=.6))
                continue
            marker = "*" if np.isfinite(fdr.iloc[i, j]) and fdr.iloc[i, j] < .05 else ("·" if np.isfinite(pvals.iloc[i, j]) and pvals.iloc[i, j] < .05 else "")
            ax.text(j+.5, i+.5, f"{value:.1f}{marker}", ha="center", va="center", fontsize=6.3,
                    color="white" if abs(value) > 1.35 else COLORS["ink"])
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_title("Human GBS module effects across independent platforms and phases", loc="left", fontweight="bold", pad=12)
    ax.tick_params(axis="x", rotation=24); ax.tick_params(axis="y", rotation=0)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.text(0, -0.34, "Cells show Hedges g; the paired acute-to-1-year plasma row shows paired dz. "
            "* within-row FDR<0.05; · nominal P<0.05. Grey = unavailable. No cross-platform meta-analysis was performed.",
            transform=ax.transAxes, fontsize=6.4, color=COLORS["gray"], va="top")
    save(fig, "Figure_2_acute_GBS_source_signatures")


def figure3() -> None:
    g211 = pd.read_csv(TABLE / "gse211225_sample_module_scores.csv")
    g211["module_label"] = g211["module"].map(MODULE_LABELS)
    selected = ["CXCL8-CXCR1/2", "LIF/OSM-gp130", "Complement", "Fc receptor", "Transmigration"]
    g211 = g211.loc[g211["module_label"].isin(selected)].copy()
    condition_label = {"Acute_GBS": "Acute GBS", "Postacute_GBS": "Post-acute GBS", "HC": "Healthy"}
    g211["stage"] = g211["condition"].map(condition_label)
    prot = pd.read_csv(TABLE / "gbs_proteomics_module_effects.csv")
    paired = prot.loc[prot["comparison"].eq("Acute_GBS_vs_Recovery_1y")].copy()
    paired["module_label"] = paired["module"].map(MODULE_LABELS)
    paired = paired.loc[paired["module_label"].isin(CORE)].sort_values("delta")
    score = pd.read_csv(TABLE / "gbs_proteomics_sample_module_scores.csv")

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), gridspec_kw={"height_ratios": [1.15, 1]}, constrained_layout=True)
    ax = axes[0]; panel(ax, "a")
    palette = {"Acute GBS": COLORS["red"], "Post-acute GBS": COLORS["orange"], "Healthy": COLORS["blue"]}
    sns.boxplot(data=g211, x="module_label", y="score_z", hue="stage", order=selected,
                palette=palette, fliersize=0, linewidth=.8, width=.72, ax=ax)
    sns.stripplot(data=g211, x="module_label", y="score_z", hue="stage", order=selected,
                  palette=palette, dodge=True, size=3.0, alpha=.75, ax=ax)
    handles, labels = ax.get_legend_handles_labels(); ax.legend(handles[:3], labels[:3], frameon=False, ncol=3, loc="upper right")
    ax.axhline(0, color=COLORS["gray"], lw=.6); ax.set_xlabel(""); ax.set_ylabel("Within-cohort module z score")
    ax.set_title("GSE211225 cross-sectional acute and post-acute whole-blood states", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=18)

    ax = axes[1]; panel(ax, "b")
    estimates=[]
    for row in paired.itertuples(index=False):
        values = score.loc[score["module"].eq(row.module)].copy()
        a = values.loc[values["condition"].eq("GBS_Acute")].set_index("matching_ID")["score_z"]
        r = values.loc[values["condition"].eq("GBS_Recovery")].set_index("matching_ID")["score_z"]
        common = a.index.intersection(r.index); diff = a.loc[common] - r.loc[common]
        sem = stats.sem(diff); ci = stats.t.ppf(.975, len(diff)-1) * sem
        estimates.append((row.module_label, float(diff.mean()), float(ci), float(row.paired_t_p), float(row.module_fdr_within_comparison)))
    est = pd.DataFrame(estimates, columns=["module", "delta", "ci", "p", "fdr"]).sort_values("delta")
    y = np.arange(len(est)); ax.errorbar(est["delta"], y, xerr=est["ci"], fmt="o", color=COLORS["purple"], ecolor=COLORS["purple"], capsize=3)
    ax.axvline(0, color=COLORS["gray"], lw=.8); ax.set_yticks(y, est["module"]); ax.set_xlabel("Acute - paired 1-year recovery module score (95% CI)")
    ax.set_title("Longitudinal plasma proteome: 20 matched GBS participants (preprint tier)", loc="left", fontweight="bold")
    for yi, row in enumerate(est.itertuples(index=False)):
        ax.text(max(est["delta"] + est["ci"]) + .06, yi, f"P={row.p:.3f}", va="center", fontsize=6.3, color=COLORS["gray"])
    save(fig, "Figure_3_GBS_phase_and_longitudinal_validation")


def figure4() -> None:
    cidp = pd.read_csv(TABLE / "gse285983_all_module_effects.csv")
    cidp = cidp.loc[cidp["comparison"].eq("CIDP_vs_CIAP")].copy()
    cidp["module_label"] = cidp["panel"].map(MODULE_LABELS)
    compartments = ["Macrophage", "BNB_EC", "Myelinating_SC", "Repair_damage_SC"]
    modules = ["CXCL8-CXCR1/2", "LIF/OSM-gp130", "Complement", "Fc receptor", "Transmigration", "BNB integrity", "Schwann/myelin repair"]
    matrix = cidp.pivot(index="cell_group", columns="module_label", values="delta_mean_z").reindex(index=compartments, columns=modules)
    fdr = cidp.pivot(index="cell_group", columns="module_label", values="fdr_within_celltype_comparison").reindex(index=compartments, columns=modules)
    bnb = pd.read_csv(TABLE / "gse107574_bnb_target_summary.csv")
    genes = ["CLDN5", "OCLN", "TJP1", "CDH5", "VWF", "LIFR", "IL6ST", "ICAM1", "VCAM1", "C3", "C5AR1"]
    bnb = bnb.loc[bnb["gene"].isin(genes)].copy(); bnb["log10_median_plus"] = np.log10(bnb["median_fpkm"] + .1)

    fig = plt.figure(figsize=(7.2, 9.5))
    fig.subplots_adjust(left=.19, right=.94, top=.96, bottom=.06, hspace=.86)
    grid = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.05, 1.05])
    ax = fig.add_subplot(grid[0]); panel(ax, "a")
    sns.heatmap(matrix, cmap="vlag", center=0, vmin=-.8, vmax=.8, linewidths=.6, linecolor="white",
                cbar_kws={"label": "CIDP - CIAP, delta mean z", "shrink": .7}, ax=ax)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix.iloc[i,j]):
                mark = "*" if fdr.iloc[i,j] < .05 else ("·" if fdr.iloc[i,j] < .10 else "")
                ax.text(j+.5, i+.5, f"{matrix.iloc[i,j]:.2f}{mark}", ha="center", va="center", fontsize=6.2,
                        color="white" if abs(matrix.iloc[i,j]) > .42 else COLORS["ink"])
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("CIDP sural nerve donor-level compartment remodeling", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=22); ax.tick_params(axis="y", rotation=0)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.text(0, -0.54, "* FDR<0.05; · FDR<0.10. CIAP is the primary neuropathy comparator.", transform=ax.transAxes, fontsize=6.3, color=COLORS["gray"])

    ax = fig.add_subplot(grid[1]); panel(ax, "b")
    prep_label = {"cultured_endoneurial_endothelial_cell": "Cultured endoneurial EC (n=2)", "LCM_endoneurial_microvessel": "LCM microvessels (n=4)"}
    bnb["Preparation"] = bnb["preparation"].map(prep_label)
    sns.scatterplot(data=bnb, x="gene", y="log10_median_plus", hue="Preparation", style="Preparation",
                    palette=[COLORS["blue"], COLORS["green"]], s=58, ax=ax)
    ax.set_xlabel(""); ax.set_ylabel("log10(median FPKM + 0.1)"); ax.tick_params(axis="x", rotation=25)
    ax.set_title("GSE107574 normal human BNB identity anchor (descriptive, not a disease contrast)", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=2, loc="upper right", fontsize=6.4)

    ax = fig.add_subplot(grid[2]); panel(ax, "c"); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 3.6)
    ax.set_title("Expression-supported myeloid-BNB-Schwann communication model", loc="left", fontweight="bold", pad=2)
    nodes = [(1.1, 1.8, "Circulating\nmonocyte", COLORS["red"]), (4.2, 1.8, "BNB\nendothelium", COLORS["blue"]),
             (7.0, 2.75, "Nerve\nmacrophage", COLORS["orange"]), (8.8, .8, "Schwann\ncell", COLORS["purple"])]
    for x,y,label,color in nodes:
        ax.add_patch(FancyBboxPatch((x-.72,y-.42),1.44,.84,boxstyle="round,pad=.05",fc="white",ec=color,lw=1.5)); ax.text(x,y,label,ha="center",va="center",fontsize=7,fontweight="bold")
    edges = [((1.85,2.05),(3.42,2.05),"CXCL8 / CCL2\nadhesion + recruitment",COLORS["red"]),
             ((4.9,2.15),(6.3,2.6),"transmigration",COLORS["green"]),
             ((7.65,2.45),(8.25,1.2),"complement / FcR\nmyelin injury-clearance",COLORS["orange"]),
             ((4.9,1.45),(8.05,.85),"LIF/OSM-gp130\nbarrier-glial adaptation",COLORS["purple"])]
    for start,end,label,color in edges:
        ax.add_patch(FancyArrowPatch(start,end,arrowstyle="-|>",mutation_scale=10,lw=1.4,color=color)); ax.text((start[0]+end[0])/2,(start[1]+end[1])/2+.18,label,ha="center",fontsize=5.9,color=color)
    ax.text(5, .05, "Edges show cross-cohort expression compatibility only; no matched donor or causal communication estimate is available.", ha="center", fontsize=6.3, color=COLORS["gray"])
    save(fig, "Figure_4_CIDP_BNB_compartment_map")


def figure5() -> None:
    ean = pd.read_csv(TABLE / "gse133750_module_effects.csv")
    ean = ean.loc[ean["comparison"].eq("Peak_neuritis_vs_Control")].copy()
    ean["module_label"] = ean["module"].map(MODULE_LABELS)
    ean = ean.loc[ean["module_label"].isin(CORE + ["Schwann/myelin repair", "BNB integrity"])].sort_values("hedges_g")
    oep = pd.read_csv(TABLE / "oep002315_oep002701_selected_author_gsea.csv")
    curated = [
        ("CD14+CD163high monocytes KEGG", "Phagosome"),
        ("CD14+CD163high monocytes KEGG", "Chemokine signaling pathway"),
        ("CD14+CD163high monocytes KEGG", "Leukocyte transendothelial migration"),
        ("CD14+MALAT1+ monocytes KEGG", "Leukocyte transendothelial migration"),
        ("CD14+MALAT1+ monocytes KEGG", "Fc gamma R-mediated phagocytosis"),
        ("CD16+ monocytes KEGG", "Complement and coagulation cascades"),
    ]
    picked=[]
    for cell, term in curated:
        row=oep.loc[oep["cell_subset"].eq(cell)&oep["description"].eq(term)]
        if len(row):
            item=row.iloc[0].copy(); item["label"]=cell.replace(" monocytes KEGG","")+" | "+term; picked.append(item)
    picked=pd.DataFrame(picked)

    phase_cols = ["Acute GBS", "Post-acute", "Paired 1-y plasma", "CIDP nerve", "Stable treated CIDP", "IVIg response"]
    phase_rows = ["CXCL8-CXCR1/2", "LIF/OSM-gp130", "Complement", "Fc receptor", "Transmigration"]
    # 0 access/absent; 1 publication-level; 2 direct descriptive; 3 direct inferential.
    evidence = np.array([
        [3,3,2,3,1,0], [3,3,3,3,1,0], [3,3,3,3,1,0], [3,3,3,3,1,1], [3,3,3,3,1,0]
    ])
    fig = plt.figure(figsize=(7.2, 9.4))
    fig.subplots_adjust(left=.32, right=.94, top=.96, bottom=.07, hspace=.78)
    grid=fig.add_gridspec(3,1,height_ratios=[1.05,1,1.15])
    ax=fig.add_subplot(grid[0]); panel(ax,"a")
    cmap=mpl.colors.ListedColormap(["#E5E7EB","#F5C995","#B9D9EA","#70AE8C"])
    sns.heatmap(evidence,cmap=cmap,vmin=0,vmax=3,cbar=False,linewidths=.7,linecolor="white",xticklabels=phase_cols,yticklabels=phase_rows,ax=ax)
    labels={0:"access gate",1:"publication",2:"descriptive",3:"direct"}
    for i in range(evidence.shape[0]):
        for j in range(evidence.shape[1]):
            ax.text(j+.5,i+.5,labels[evidence[i,j]],ha="center",va="center",fontsize=5.8,color="white" if evidence[i,j]==3 else COLORS["ink"])
    ax.tick_params(axis="x",rotation=22); ax.tick_params(axis="y",rotation=0); ax.set_xlabel("");ax.set_ylabel("")
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.set_title("Phase and treatment coverage after the upgrade",loc="left",fontweight="bold")

    ax=fig.add_subplot(grid[1]); panel(ax,"b")
    ax.barh(ean["module_label"],ean["hedges_g"],color=[COLORS["red"] if v>0 else COLORS["blue"] for v in ean["hedges_g"]]);ax.axvline(0,color=COLORS["gray"],lw=.7)
    ax.set_xlabel("Rat EAN peak neuritis vs control, Hedges g");ax.set_title("GSE133750 sciatic-nerve EAN mechanism tier (n=3/group)",loc="left",fontweight="bold")

    ax=fig.add_subplot(grid[2]); panel(ax,"c")
    picked["short_label"] = picked["label"].replace({
        "CD16+ | Complement and coagulation cascades": "CD16+ | complement/coagulation",
        "CD14+MALAT1+ | Fc gamma R-mediated phagocytosis": "CD14+MALAT1+ | FcγR phagocytosis",
        "CD14+MALAT1+ | Leukocyte transendothelial migration": "CD14+MALAT1+ | transmigration",
        "CD14+CD163high | Leukocyte transendothelial migration": "CD14+CD163hi | transmigration",
        "CD14+CD163high | Chemokine signaling pathway": "CD14+CD163hi | chemokine signaling",
        "CD14+CD163high | Phagosome": "CD14+CD163hi | phagosome",
    })
    ax.barh(picked["short_label"],picked["NES"],color=COLORS["orange"]);ax.axvline(0,color=COLORS["gray"],lw=.7)
    ax.set_xlabel("Author-reported normalized enrichment score");ax.set_title("OEP002315/OEP002701 monocyte external validation (supplementary GSEA)",loc="left",fontweight="bold")
    for yi,row in enumerate(picked.itertuples(index=False)):
        ax.text(row.NES+.04,yi,f"FDR={row.p_adjust:.3f}",va="center",fontsize=6.1,color=COLORS["gray"])
    save(fig,"Figure_5_phase_aware_continuum")


def main() -> None:
    style()
    effects = module_effect_matrix()
    stats_payload = key_stats(effects)
    figure1(); figure2(effects.copy()); figure3(); figure4(); figure5()
    print(f"Wrote {OUT_JSON}")
    print(f"Direct human resources: {stats_payload['direct_human_resource_count']}")
    print(f"Wrote figures to {FIG}")


if __name__ == "__main__":
    main()
