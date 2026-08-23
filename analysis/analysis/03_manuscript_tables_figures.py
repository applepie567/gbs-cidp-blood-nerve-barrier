#!/usr/bin/env python3
"""Build manuscript-ready tables and figures for the GBS-CIDP BNB study."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results/tables"
FIG_DIR = ROOT / "submission/figures"
CIDP_JSON = ROOT / "results/gse285983_targeted_pseudobulk.json"
GBS_GENES = TABLE_DIR / "gse304871_targeted_gene_effects.csv"
GBS_MODULES = TABLE_DIR / "gse304871_module_effects.csv"
GBS_SCORES = TABLE_DIR / "gse304871_sample_module_scores.csv"
CIDP_EXPR = TABLE_DIR / "gse285983_sample_targeted_expression.csv"
CIDP_MODULES = TABLE_DIR / "gse285983_sample_module_scores.csv"


COLORS = {
    "blue": "#2878B5",
    "orange": "#E07A2D",
    "green": "#3A9D72",
    "purple": "#7A5AA6",
    "red": "#C94C4C",
    "gray": "#6B7280",
    "light_blue": "#DCEAF5",
    "light_orange": "#F8E3D3",
    "light_green": "#DDEFE7",
    "light_purple": "#E9E2F2",
    "ink": "#20242A",
}


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    sns.set_style("whitegrid", {"grid.color": "#E5E7EB", "grid.linewidth": 0.5})


def panel_label(axis, label: str) -> None:
    axis.text(-0.08, 1.04, label, transform=axis.transAxes, fontsize=11,
              fontweight="bold", va="bottom", ha="left", color=COLORS["ink"])


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def build_tables(cidp: dict, gbs_genes: pd.DataFrame, gbs_modules: pd.DataFrame) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    datasets = pd.DataFrame([
        {
            "role": "Primary acute immune-source dataset",
            "accession": "GSE304871",
            "disease_stage": "Early untreated AIDP-variant GBS",
            "biospecimen": "Peripheral blood; FACS-sorted CD11b+, CD4+, CD8+ cells",
            "biological_replicates": "CD11b: 2 GBS/3 HC; CD4: 3/3; CD8: 2/3",
            "analysis": "Patient-level log2(nTPM+1), prespecified modules and Hedges g",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304871",
        },
        {
            "role": "Primary target-organ dataset",
            "accession": "GSE285983",
            "disease_stage": "Cross-sectional CIDP nerve biopsy",
            "biospecimen": "Human sural nerve single-nucleus RNA-seq",
            "biological_replicates": "9 CIDP; 11 CIAP; 4 trauma-graft controls; 37 total donors",
            "analysis": "Patient-level targeted pseudobulk; CIDP-vs-CIAP primary contrast",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285983",
        },
        {
            "role": "External phase/treatment validation",
            "accession": "Published cohorts",
            "disease_stage": "GBS pre/post treatment; CIDP active/inactive and IVIg response",
            "biospecimen": "CSF, serum, whole blood and cultured PnMECs",
            "biological_replicates": "Study-specific",
            "analysis": "Directional triangulation only; no unreported individual-level data pooled",
            "source_url": "https://doi.org/10.3389/fimmu.2023.1241199",
        },
    ])
    datasets.to_csv(TABLE_DIR / "table1_datasets_and_stage_coverage.csv", index=False)

    patient_info = pd.read_excel(
        ROOT / "data/raw/GSE285983/supplement/MOESM3.xlsx",
        sheet_name="patient_info",
    )
    patient_info.columns = [str(column).strip().replace(" ", "_").replace("µ", "u") for column in patient_info.columns]
    patient_info.to_csv(TABLE_DIR / "gse285983_patient_metadata.csv", index=False)

    comparisons = pd.DataFrame(cidp["comparisons"])
    atlas = pd.DataFrame(cidp["expression_atlas"])
    pooled_atlas = pd.DataFrame(cidp["pooled_expression_atlas"])
    module = pd.DataFrame(cidp["module_comparisons"])
    abundance = pd.DataFrame(cidp["abundance_comparisons"])
    comparisons.to_csv(TABLE_DIR / "gse285983_all_targeted_gene_effects.csv", index=False)
    atlas.to_csv(TABLE_DIR / "gse285983_expression_atlas.csv", index=False)
    pooled_atlas.to_csv(TABLE_DIR / "gse285983_pooled_expression_atlas.csv", index=False)
    module.to_csv(TABLE_DIR / "gse285983_all_module_effects.csv", index=False)
    abundance.to_csv(TABLE_DIR / "gse285983_cell_abundance_effects.csv", index=False)

    direct_rows = []
    axes = {
        "CXCL8-CXCR1/2": ("CXCL8", "Macrophage", "CXCR1;CXCR2", "Granulocyte"),
        "LIF-LIFR-gp130": ("LIF", "Macrophage", "LIFR;IL6ST", "BNB_EC"),
        "OSM-OSMR-gp130": ("OSM", "Macrophage", "OSMR;IL6ST", "BNB_EC"),
        "CCL2-CCR2": ("CCL2", "Macrophage", "CCR2", "Macrophage"),
        "CCL20-CCR6": ("CCL20", "BNB_EC", "CCR6", "T_NK"),
        "C3-C3AR1": ("C3", "Macrophage", "C3AR1", "Macrophage"),
        "C5-C5AR1": ("C5", "Macrophage", "C5AR1", "Macrophage"),
        "SELP-SELPLG": ("SELP", "BNB_EC", "SELPLG", "Macrophage"),
    }
    for axis, (ligand, sender, receptors, receiver) in axes.items():
        g = gbs_genes[(gbs_genes.cell_type == "CD11b") & (gbs_genes.gene == ligand)]
        ligand_delta = float(g.delta_log2.iloc[0]) if len(g) else np.nan
        ligand_g = float(g.hedges_g.iloc[0]) if len(g) else np.nan
        receptor_values = []
        for receptor in receptors.split(";"):
            row = pooled_atlas[(pooled_atlas.disease == "CIDP") & (pooled_atlas.cell_group == receiver) & (pooled_atlas.gene == receptor)]
            if len(row):
                receptor_values.append(float(row.pooled_log2cpm.iloc[0]))
        direct_rows.append({
            "axis": axis,
            "acute_GBS_CD11b_ligand": ligand,
            "acute_GBS_delta_log2_nTPM": ligand_delta,
            "acute_GBS_hedges_g": ligand_g,
            "CIDP_nerve_sender": sender,
            "CIDP_nerve_receptors": receptors,
            "CIDP_nerve_receiver": receiver,
            "CIDP_mean_receptor_log2CPM": float(np.mean(receptor_values)) if receptor_values else np.nan,
            "interpretation": "Expression-supported, hypothesis-generating; not a causal interaction estimate",
        })
    pd.DataFrame(direct_rows).to_csv(TABLE_DIR / "table2_cross_compartment_axes.csv", index=False)

    phase = pd.DataFrame([
        ["CXCL8-CXCR1/2", "Direct: CXCL8 higher in acute GBS CD11b+ cells", "Published: CSF IL8 highest in pretreatment GBS", "Published: CSF IL8 falls after GBS treatment", "Direct: macrophage-localized CXCL8; granulocyte CXCR1/2", "Published: CIDP CSF IL8 remains high after treatment", "https://doi.org/10.3389/fimmu.2023.1241199"],
        ["LIF-LIFR-gp130", "Direct: OSM stronger than LIF in acute GBS CD11b+ cells", "Published: CSF LIF elevated in GBS", "No public recovery matrix", "Direct: LIFR/IL6ST high in BNB endothelium and Schwann cells", "No paired nerve data", "https://doi.org/10.1007/s00109-026-02698-2"],
        ["Complement", "Direct: C3 and CFB directionally higher in acute GBS CD11b+ cells", "No harmonized phase estimate", "Published paired serum proteomics supports acute-to-recovery remodeling", "Direct: macrophage C3/C3AR1 and complement-regulatory expression", "IVIg can modulate complement; phase-specific raw data unavailable", "https://doi.org/10.64898/2026.05.23.26353948"],
        ["Fc receptor", "Direct: mixed acute GBS CD11b+ Fc-receptor changes", "Not measured in CSF validation", "Not measured", "Direct: macrophage-dominant Fc genes; module lower vs CIAP", "Published: IVIg partially restores Fc-gamma receptor balance", "https://doi.org/10.1212/NXI.0000000000000148"],
        ["Leukocyte transmigration", "Direct: CCR1/CCL2/CCL20 directionally higher in acute GBS CD11b+ cells", "Published: CSF SELE/ITGAM signals in GBS", "No public paired transcriptome", "Direct: BNB ICAM1/SELP/VCAM1 competence", "Published: CIDP IgG increases BNB permeability and endothelial CCL20/VCAM1", "https://doi.org/10.3390/ijms27021088"],
        ["CIDP activity", "Not applicable", "Not applicable", "Not applicable", "Cross-sectional tissue only", "Published: IRAK4 plus four proteins tracked active/inactive change", "https://doi.org/10.1136/jnnp-2023-332398"],
        ["IVIg response", "Not applicable", "Not applicable", "Not applicable", "Biopsy metadata include response label, not paired dosing", "Published: TNFR1/TLR and Fc-receptor programs change after IVIg", "https://doi.org/10.1097/MD.0000000000003370"],
    ], columns=["axis", "GBS_acute_blood", "GBS_acute_CSF", "GBS_recovery_or_post_treatment", "CIDP_nerve", "CIDP_active_stable_or_IVIg", "source_url"])
    phase.to_csv(TABLE_DIR / "phase_evidence_matrix.csv", index=False)

    summary = {
        "cidp_n": sum(1 for value in cidp["sample_groups"].values() if value == "CIDP"),
        "ciap_n": sum(1 for value in cidp["sample_groups"].values() if value == "CIAP"),
        "ctrl_n": sum(1 for value in cidp["sample_groups"].values() if value == "CTRL"),
        "total_nuclei": int(sum(cidp["total_nuclei"].values())),
        "gbs_cd11b_n": "2 GBS vs 3 HC",
        "cidp_macrophage_fc_module": module[(module.comparison == "CIDP_vs_CIAP") & (module.cell_group == "Macrophage") & (module.panel == "Fc_receptor")].to_dict("records")[0],
        "abundance_min_fdr": float(abundance.fdr_across_cell_groups.min()),
    }
    (ROOT / "results/manuscript_stats.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def figure1() -> None:
    fig = plt.figure(figsize=(7.1, 6.1))
    grid = fig.add_gridspec(2, 1, height_ratios=[1.65, 1], hspace=0.28)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    panel_label(ax, "a")
    ax.text(0, 6.15, "Compartment-aware myeloid-vascular-glial model", fontsize=10, fontweight="bold")

    ax.add_patch(FancyBboxPatch((0.15, 0.45), 2.25, 4.7, boxstyle="round,pad=0.08",
                                facecolor=COLORS["light_blue"], edgecolor=COLORS["blue"], lw=1.2))
    ax.text(1.28, 4.82, "Blood\nacute GBS", ha="center", va="top", fontweight="bold", color=COLORS["blue"])
    for y, label, color in [(3.55, "CD11b+ myeloid", COLORS["orange"]), (2.55, "CD4+ T cells", COLORS["green"]), (1.55, "CD8+ T cells", COLORS["purple"])]:
        ax.add_patch(plt.Circle((1.28, y), 0.36, fc=color, ec="white", lw=1.2))
        ax.text(1.28, y-0.62, label, ha="center", va="top", fontsize=7)

    ax.add_patch(FancyBboxPatch((3.0, 0.45), 2.65, 4.7, boxstyle="round,pad=0.08",
                                facecolor="#EEF3F7", edgecolor="#4D6A7F", lw=1.2))
    ax.text(4.33, 4.82, "Blood-nerve barrier (BNB)", ha="center", va="top", fontweight="bold")
    for x in [3.45, 4.05, 4.65, 5.25]:
        ax.add_patch(Rectangle((x-0.25, 2.45), 0.5, 0.72, fc=COLORS["blue"], ec="white", lw=0.8))
    ax.text(4.33, 2.32, "endoneurial microvascular\nendothelium + pericytes", ha="center", va="top", fontsize=6.6)
    ax.text(4.33, 0.72, "BNB is distinct from the BBB", ha="center", va="center", color=COLORS["red"], fontweight="bold")

    ax.add_patch(FancyBboxPatch((6.25, 0.45), 3.55, 4.7, boxstyle="round,pad=0.08",
                                facecolor=COLORS["light_green"], edgecolor=COLORS["green"], lw=1.2))
    ax.text(8.02, 4.82, "Peripheral nerve endoneurium\nCIDP target-organ state", ha="center", va="top", fontweight="bold")
    ax.add_patch(plt.Circle((7.25, 2.85), 0.48, fc=COLORS["orange"], ec="white"))
    ax.text(7.25, 2.1, "macrophage", ha="center", fontsize=7)
    ax.add_patch(FancyBboxPatch((8.25, 2.4), 1.05, 0.9, boxstyle="round,pad=0.04",
                                fc=COLORS["purple"], ec="white"))
    ax.text(8.78, 2.02, "Schwann cell", ha="center", fontsize=7)
    ax.plot([8.78, 8.78], [3.3, 4.0], color=COLORS["purple"], lw=5, solid_capstyle="round")
    ax.text(8.78, 1.1, "myelin / repair", ha="center", fontsize=7)

    for y, text_label, color in [
        (4.05, "CXCL8 -> CXCR1/2", COLORS["red"]),
        (3.55, "LIF/OSM -> LIFR/OSMR-gp130", COLORS["purple"]),
        (1.50, "complement + Fc receptors", COLORS["orange"]),
        (1.05, "adhesion and transmigration", COLORS["blue"]),
    ]:
        ax.add_patch(FancyArrowPatch((2.42, y), (6.2, y), arrowstyle="-|>", mutation_scale=10,
                                     lw=1.5, color=color))
        ax.text(4.32, y+0.10, text_label, ha="center", va="bottom", fontsize=6.3, color=color)

    ax2 = fig.add_subplot(grid[1, 0])
    ax2.axis("off")
    panel_label(ax2, "b")
    ax2.text(0, 1.03, "Phase-aware evidence architecture", transform=ax2.transAxes, fontsize=10, fontweight="bold")
    stages = [
        (0.03, "GBS\nacute", "Direct reanalysis\nGSE304871", COLORS["red"]),
        (0.27, "GBS post-treatment\n/ recovery", "Published longitudinal\nCSF/proteomics", COLORS["orange"]),
        (0.56, "CIDP\nnerve state", "Direct reanalysis\nGSE285983", COLORS["blue"]),
        (0.80, "CIDP active / stable\n/ IVIg", "Published longitudinal\nserum/blood", COLORS["green"]),
    ]
    for x, title, subtitle, color in stages:
        ax2.add_patch(FancyBboxPatch((x, 0.25), 0.18, 0.5, boxstyle="round,pad=0.02",
                                     transform=ax2.transAxes, fc="white", ec=color, lw=1.4))
        ax2.text(x+0.09, 0.61, title, transform=ax2.transAxes, ha="center", va="center", fontweight="bold", fontsize=6.4)
        ax2.text(x+0.09, 0.39, subtitle, transform=ax2.transAxes, ha="center", va="center", fontsize=6.5)
    for x in [0.215, 0.505, 0.745]:
        ax2.annotate("", xy=(x+0.045, 0.5), xytext=(x, 0.5), xycoords=ax2.transAxes,
                     arrowprops=dict(arrowstyle="->", color=COLORS["gray"], lw=1.1))
    ax2.text(0.5, 0.05, "Biological samples, not individual cells/nuclei, are the inferential units",
             transform=ax2.transAxes, ha="center", fontsize=7, color=COLORS["gray"])
    save_figure(fig, "Figure_1_study_design")


def figure2(gbs_genes: pd.DataFrame, gbs_scores: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 7.0), gridspec_kw={"height_ratios": [1.15, 1]}, constrained_layout=True)
    ax = axes[0]
    panel_label(ax, "a")
    selection = ["OSM", "OAS1", "CD59", "CCL2", "CXCL8", "CCR1", "C3", "IFI6", "IL6R", "CFB", "FCGR2A", "SOCS3"]
    data = gbs_genes[(gbs_genes.cell_type == "CD11b") & (gbs_genes.gene.isin(selection))].copy()
    data = data.sort_values("delta_log2")
    axis_map = {
        "OSM": "IL-6 family", "IL6R": "IL-6 family", "SOCS3": "IL-6 family",
        "CXCL8": "CXCL8", "CCR1": "Migration", "CCL2": "Migration",
        "C3": "Complement", "CFB": "Complement", "CD59": "Complement",
        "FCGR2A": "Fc receptor", "IFI6": "Interferon", "OAS1": "Interferon",
    }
    palette = {"IL-6 family": COLORS["purple"], "CXCL8": COLORS["red"], "Migration": COLORS["blue"],
               "Complement": COLORS["orange"], "Fc receptor": COLORS["green"], "Interferon": COLORS["gray"]}
    bar_colors = [palette[axis_map[g]] for g in data.gene]
    ax.barh(data.gene, data.delta_log2, color=bar_colors, edgecolor="white")
    ax.axvline(0, color=COLORS["ink"], lw=0.8)
    ax.set_xlabel("GBS - healthy control, mean log2(nTPM + 1)")
    ax.set_ylabel("")
    ax.set_title("Early untreated GBS CD11b+ cells (2 GBS vs 3 controls)", loc="left", fontweight="bold")
    ax.text(0.99, 0.02, "Effect-size screen; no targeted FDR < 0.05", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.5, color=COLORS["gray"])

    ax = axes[1]
    panel_label(ax, "b")
    modules = ["CXCL8_axis", "IL6_family_LIF_axis", "Transendothelial_migration", "Complement", "Fc_receptor"]
    scores = gbs_scores[(gbs_scores.cell_type == "CD11b") & (gbs_scores.module.isin(modules))].copy()
    label_map = {"CXCL8_axis": "CXCL8", "IL6_family_LIF_axis": "LIF/OSM-gp130",
                 "Transendothelial_migration": "Migration", "Complement": "Complement", "Fc_receptor": "Fc receptor"}
    scores["module_label"] = scores.module.map(label_map)
    order = list(label_map.values())
    sns.boxplot(data=scores, x="module_label", y="score_z", hue="group", order=order,
                palette={"GBS": COLORS["red"], "HC": COLORS["light_blue"]}, width=0.55,
                fliersize=0, linewidth=0.8, ax=ax)
    sns.stripplot(data=scores, x="module_label", y="score_z", hue="group", order=order,
                  dodge=True, palette={"GBS": COLORS["red"], "HC": COLORS["blue"]},
                  size=4, linewidth=0.4, edgecolor="white", ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="", frameon=False, loc="upper right")
    ax.axhline(0, color=COLORS["gray"], lw=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Within-fraction module z score")
    ax.set_title("Patient-level prespecified module scores", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=15)
    save_figure(fig, "Figure_2_acute_GBS_source_signatures")


def figure3(cidp: dict, cidp_expr: pd.DataFrame) -> None:
    atlas = pd.DataFrame(cidp["pooled_expression_atlas"])
    fig = plt.figure(figsize=(7.1, 8.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.1, 1])
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "a")
    cell_groups = ["Macrophage", "Granulocyte", "BNB_EC", "Other_EC", "Pericyte", "Perineurium", "Myelinating_SC", "Repair_damage_SC"]
    genes = ["CXCL8", "CXCR1", "CXCR2", "LIFR", "IL6ST", "OSMR", "C3", "C3AR1", "FCGR2A", "FCGRT",
             "ICAM1", "VCAM1", "SELP", "ABCB1", "MFSD2A", "CLDN5", "MPZ", "NGFR"]
    data = atlas[(atlas.disease == "CIDP") & (atlas.cell_group.isin(cell_groups)) & (atlas.gene.isin(genes))]
    matrix = data.pivot(index="gene", columns="cell_group", values="pooled_log2cpm").reindex(index=genes, columns=cell_groups)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0)
    sns.heatmap(z, cmap="vlag", center=0, vmin=-2, vmax=2, ax=ax, cbar_kws={"label": "Gene-wise z score", "shrink": 0.65},
                linewidths=0.35, linecolor="white")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("CIDP nerve: compartment localization of prespecified axes", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=32)
    ax.tick_params(axis="y", rotation=0)

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "b")
    landmarks = [
        ("Macrophage", "CXCL8", "Macrophage CXCL8"),
        ("BNB_EC", "LIFR", "BNB LIFR"),
        ("BNB_EC", "ABCB1", "BNB ABCB1"),
        ("Macrophage", "C3", "Macrophage C3"),
        ("Macrophage", "FCGR2A", "Macrophage FCGR2A"),
        ("Repair_damage_SC", "NGFR", "Repair/damage SC NGFR"),
    ]
    pieces = []
    for cell_group, gene, label in landmarks:
        subset = cidp_expr[(cidp_expr.cell_group == cell_group) & (cidp_expr.gene == gene) & (cidp_expr.disease.isin(["CIDP", "CIAP"]))].copy()
        subset["landmark"] = label
        pieces.append(subset)
    plot = pd.concat(pieces, ignore_index=True)
    order = [item[2] for item in landmarks]
    sns.boxplot(data=plot, x="landmark", y="log2_cpm_plus_0_5", hue="disease", order=order,
                palette={"CIDP": COLORS["blue"], "CIAP": COLORS["light_orange"]}, width=0.6,
                fliersize=0, linewidth=0.8, ax=ax)
    sns.stripplot(data=plot, x="landmark", y="log2_cpm_plus_0_5", hue="disease", order=order,
                  dodge=True, palette={"CIDP": COLORS["blue"], "CIAP": COLORS["orange"]},
                  size=3.3, linewidth=0.3, edgecolor="white", ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="", frameon=False, loc="upper right")
    ax.set_xlabel("")
    ax.set_ylabel("Patient pseudobulk log2(CPM + 0.5)")
    ax.set_title("Expression competence is stronger than disease-specific induction", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=24)
    save_figure(fig, "Figure_3_CIDP_BNB_compartment_map")


def figure4(cidp: dict, gbs_genes: pd.DataFrame) -> None:
    atlas = pd.DataFrame(cidp["pooled_expression_atlas"])
    fig = plt.figure(figsize=(7.1, 7.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 1.15])
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "a")
    pairs = [
        ("CXCL8", "CXCR1", "Granulocyte", "CXCL8-CXCR1"),
        ("CXCL8", "CXCR2", "Granulocyte", "CXCL8-CXCR2"),
        ("OSM", "OSMR", "BNB_EC", "OSM-OSMR"),
        ("LIF", "LIFR", "BNB_EC", "LIF-LIFR"),
        ("CCL2", "CCR2", "Macrophage", "CCL2-CCR2"),
        ("CCL20", "CCR6", "T_NK", "CCL20-CCR6"),
        ("C3", "C3AR1", "Macrophage", "C3-C3AR1"),
        ("C5", "C5AR1", "Macrophage", "C5-C5AR1"),
    ]
    rows = []
    for ligand, receptor, receiver, label in pairs:
        g = gbs_genes[(gbs_genes.cell_type == "CD11b") & (gbs_genes.gene == ligand)]
        n = atlas[(atlas.disease == "CIDP") & (atlas.cell_group == receiver) & (atlas.gene == receptor)]
        if len(g) and len(n):
            rows.append({"label": label, "ligand_delta": float(g.delta_log2.iloc[0]),
                         "receptor_expression": float(n.pooled_log2cpm.iloc[0]), "receiver": receiver})
    data = pd.DataFrame(rows)
    receiver_colors = {"Granulocyte": COLORS["red"], "BNB_EC": COLORS["blue"],
                       "Macrophage": COLORS["orange"], "T_NK": COLORS["green"]}
    for receiver, subset in data.groupby("receiver"):
        ax.scatter(subset.ligand_delta, subset.receptor_expression, s=48,
                   color=receiver_colors[receiver], edgecolor="white", linewidth=0.7, label=receiver)
        for _, row in subset.iterrows():
            ax.text(row.ligand_delta + 0.035, row.receptor_expression + 0.10, row.label, fontsize=6.3)
    ax.axvline(0, color=COLORS["gray"], lw=0.7, ls="--")
    ax.set_xlabel("Acute GBS CD11b+ ligand change, log2(nTPM + 1)")
    ax.set_ylabel("CIDP nerve receiver expression, mean log2(CPM + 0.5)")
    ax.set_title("Cross-cohort ligand-receptor concordance", loc="left", fontweight="bold")
    ax.legend(title="CIDP receiver", frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.20))

    ax = fig.add_subplot(grid[1, 0])
    ax.axis("off")
    panel_label(ax, "b")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.text(0, 5.15, "Expression-supported communication model", fontsize=9, fontweight="bold")
    nodes = {
        "Circulating\nmyeloid": (1.2, 2.5, COLORS["red"]),
        "BNB\nendothelium": (5.0, 3.7, COLORS["blue"]),
        "Nerve\nmacrophage": (5.0, 1.25, COLORS["orange"]),
        "Schwann\ncell": (8.8, 2.5, COLORS["purple"]),
    }
    for label, (x, y, color) in nodes.items():
        ax.add_patch(FancyBboxPatch((x-0.75, y-0.45), 1.5, 0.9, boxstyle="round,pad=0.08",
                                    fc="white", ec=color, lw=1.6))
        ax.text(x, y, label, ha="center", va="center", fontweight="bold", fontsize=7)
    edges = [
        ((1.95, 2.8), (4.18, 3.55), "CCL2/CCL20\nadhesion competence", COLORS["blue"]),
        ((1.95, 2.2), (4.18, 1.45), "CXCL8/OSM\nmyeloid recruitment", COLORS["red"]),
        ((5.75, 3.55), (8.0, 2.8), "LIFR-gp130\nvascular-to-glial response", COLORS["purple"]),
        ((5.75, 1.45), (8.0, 2.2), "complement/Fc\nmyelin injury-repair", COLORS["orange"]),
        ((5.0, 3.15), (5.0, 1.75), "transendothelial\nmigration", COLORS["green"]),
    ]
    for start, end, label, color in edges:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10, lw=1.5, color=color))
        mx, my = (start[0]+end[0])/2, (start[1]+end[1])/2
        ax.text(mx, my+0.12, label, ha="center", va="bottom", fontsize=6.2, color=color)
    ax.text(5, 0.15, "Edges denote compatible expression and external functional evidence, not inferred causality.",
            ha="center", fontsize=6.7, color=COLORS["gray"])
    save_figure(fig, "Figure_4_cross_compartment_communication")


def figure5() -> None:
    fig, ax = plt.subplots(figsize=(7.1, 5.3))
    panel_label(ax, "a")
    rows = ["CXCL8-CXCR1/2", "LIF/OSM-gp130", "Complement", "Fc receptor", "Transmigration", "Activity markers"]
    cols = ["GBS acute\nblood", "GBS acute\nCSF", "GBS post-treatment\n/recovery", "CIDP\nnerve", "CIDP active\nvs stable", "CIDP after\nIVIg"]
    # 0 absent/not evaluated; 1 external; 2 direct; 3 direct plus orthogonal external support.
    matrix = np.array([
        [3, 3, 1, 3, 1, 1],
        [2, 3, 0, 3, 0, 0],
        [2, 0, 1, 3, 0, 1],
        [2, 0, 0, 3, 0, 3],
        [3, 1, 0, 3, 0, 2],
        [0, 0, 0, 0, 3, 2],
    ])
    cmap = mpl.colors.ListedColormap(["#F3F4F6", "#FAD9BF", "#BFDCEC", "#74B59A"])
    sns.heatmap(matrix, cmap=cmap, vmin=0, vmax=3, cbar=False, linewidths=1, linecolor="white",
                xticklabels=cols, yticklabels=rows, square=False, ax=ax)
    ax.set_title("Phase-aware evidence matrix", loc="left", fontweight="bold", pad=12)
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)
    labels = {0: "not evaluated", 1: "published", 2: "direct", 3: "triangulated"}
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            ax.text(j+0.5, i+0.5, labels[value], ha="center", va="center", fontsize=5.7,
                    color="white" if value == 3 else COLORS["ink"])
    ax.text(0, -0.18,
            "Direct = reanalysed public patient-level data; published = phase information from independent cohorts;\n"
            "triangulated = direct cross-compartment signal plus orthogonal protein/functional support.",
            transform=ax.transAxes, fontsize=6.5, va="top", color=COLORS["gray"])
    save_figure(fig, "Figure_5_phase_aware_continuum")


def main() -> None:
    set_style()
    cidp = json.loads(CIDP_JSON.read_text(encoding="utf-8"))
    gbs_genes = pd.read_csv(GBS_GENES)
    gbs_modules = pd.read_csv(GBS_MODULES)
    gbs_scores = pd.read_csv(GBS_SCORES)
    cidp_expr = pd.read_csv(CIDP_EXPR)
    build_tables(cidp, gbs_genes, gbs_modules)
    figure1()
    figure2(gbs_genes, gbs_scores)
    figure3(cidp, cidp_expr)
    figure4(cidp, gbs_genes)
    figure5()
    print(f"Wrote tables to {TABLE_DIR}")
    print(f"Wrote figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
