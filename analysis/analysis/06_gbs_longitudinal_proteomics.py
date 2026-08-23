#!/usr/bin/env python3
"""Targeted reanalysis of the public longitudinal GBS plasma proteome.

The SomaScan cohort comprises 20 matched acute/recovery pairs and 15 healthy
controls. Acute-versus-recovery tests are paired by matching_ID; control
contrasts are sample-level unpaired tests. The source is a public preprint and
is therefore retained as a clearly labelled evidence tier.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from common import GENE_MODULES, bh_adjust, paired_effect, unpaired_effect  # noqa: E402


DATA_ROOT = ROOT / "data/raw/GBS-Proteomics/Data"
TABLE_ROOT = ROOT / "results/tables"
OUT_JSON = ROOT / "results/gbs_longitudinal_proteomics.json"
ANALYSIS_MODULES = {
    key: value
    for key, value in GENE_MODULES.items()
    if key not in {"BNB_integrity", "Schwann_myelin_repair"}
}


def main() -> None:
    raw = pd.read_csv(DATA_ROOT / "GBS_Proteomics_Raw_Github.csv")
    dictionary = pd.read_csv(DATA_ROOT / "SomaScan_Protein_Dictionary_Github.csv")
    aptamers = [column for column in raw.columns if column.startswith("seq.")]
    log_rfu = np.log2(raw.set_index("SampleId")[aptamers].apply(pd.to_numeric, errors="coerce"))
    design = raw.set_index("SampleId")[["Group", "matching_ID", "Age", "Sex"]].copy()

    dictionary = dictionary.loc[dictionary["AptName"].isin(aptamers)].copy()
    dictionary["gene"] = dictionary["EntrezGeneSymbol"].fillna("").astype(str).str.split("[;, ]+")
    dictionary = dictionary.explode("gene")
    dictionary["gene"] = dictionary["gene"].str.strip()
    target_genes = sorted({gene for genes in ANALYSIS_MODULES.values() for gene in genes})
    target_map = dictionary.loc[dictionary["gene"].isin(target_genes), ["AptName", "gene", "Target", "UniProt"]]

    acute = design.index[design["Group"].eq("GBS_Acute")].tolist()
    recovery = design.index[design["Group"].eq("GBS_Recovery")].tolist()
    controls = design.index[design["Group"].eq("Healthy_Control")].tolist()
    acute_by_pair = design.loc[acute].reset_index().set_index("matching_ID")["SampleId"]
    recovery_by_pair = design.loc[recovery].reset_index().set_index("matching_ID")["SampleId"]
    pair_ids = sorted(set(acute_by_pair.index).intersection(recovery_by_pair.index))
    if len(pair_ids) != 20:
        raise ValueError(f"Expected 20 acute/recovery pairs, found {len(pair_ids)}")
    paired_acute = [int(acute_by_pair[pair]) for pair in pair_ids]
    paired_recovery = [int(recovery_by_pair[pair]) for pair in pair_ids]

    aptamer_rows: list[dict] = []
    for row in target_map.itertuples(index=False):
        values = log_rfu[row.AptName]
        tests = {
            "Acute_GBS_vs_Recovery_1y": paired_effect(values.loc[paired_acute], values.loc[paired_recovery]),
            "Acute_GBS_vs_HC": unpaired_effect(values.loc[acute], values.loc[controls]),
            "Recovery_1y_vs_HC": unpaired_effect(values.loc[recovery], values.loc[controls]),
        }
        for comparison, effect in tests.items():
            aptamer_rows.append({
                "dataset": "GBS-Proteomics",
                "comparison": comparison,
                "aptamer": row.AptName,
                "gene": row.gene,
                "target": row.Target,
                "uniprot": row.UniProt,
                "modules": ";".join(name for name, genes in ANALYSIS_MODULES.items() if row.gene in genes),
                **effect,
            })
    aptamer_effects = pd.DataFrame(aptamer_rows)
    aptamer_effects["primary_p"] = np.where(
        aptamer_effects["comparison"].eq("Acute_GBS_vs_Recovery_1y"),
        aptamer_effects["paired_t_p"],
        aptamer_effects["welch_p"],
    )
    for comparison, index in aptamer_effects.groupby("comparison").groups.items():
        aptamer_effects.loc[index, "targeted_fdr_within_comparison"] = bh_adjust(
            aptamer_effects.loc[index, "primary_p"]
        )

    # Standardize each aptamer over all 55 samples, average aptamers mapping to
    # the same gene, then average genes within each prespecified module.
    target_aptamers = target_map["AptName"].drop_duplicates().tolist()
    target_values = log_rfu[target_aptamers]
    z = target_values.sub(target_values.mean(axis=0), axis=1).div(target_values.std(axis=0, ddof=1), axis=1)
    gene_scores = {}
    for gene, subset in target_map.groupby("gene"):
        gene_scores[gene] = z[subset["AptName"].drop_duplicates()].mean(axis=1)
    gene_scores = pd.DataFrame(gene_scores)
    module_scores = {}
    availability = {}
    for module, genes in ANALYSIS_MODULES.items():
        present = [gene for gene in genes if gene in gene_scores.columns]
        availability[module] = present
        # The three-member chemokine axis can be represented by a single
        # measured circulating protein on SomaScan; n_genes is retained so
        # that this is never mistaken for a multi-protein pathway score.
        minimum = 1 if module == "CXCL8_CXCR1_2" else 2
        if len(present) >= minimum:
            module_scores[module] = gene_scores[present].mean(axis=1)
    module_scores = pd.DataFrame(module_scores)

    module_rows: list[dict] = []
    for module in module_scores.columns:
        values = module_scores[module]
        comparisons = {
            "Acute_GBS_vs_Recovery_1y": paired_effect(values.loc[paired_acute], values.loc[paired_recovery]),
            "Acute_GBS_vs_HC": unpaired_effect(values.loc[acute], values.loc[controls]),
            "Recovery_1y_vs_HC": unpaired_effect(values.loc[recovery], values.loc[controls]),
        }
        for comparison, effect in comparisons.items():
            module_rows.append({
                "dataset": "GBS-Proteomics",
                "comparison": comparison,
                "module": module,
                "n_genes": len(availability[module]),
                "genes_available": ";".join(availability[module]),
                **effect,
            })
    module_effects = pd.DataFrame(module_rows)
    module_effects["primary_p"] = np.where(
        module_effects["comparison"].eq("Acute_GBS_vs_Recovery_1y"),
        module_effects["paired_t_p"],
        module_effects["welch_p"],
    )
    for comparison, index in module_effects.groupby("comparison").groups.items():
        module_effects.loc[index, "module_fdr_within_comparison"] = bh_adjust(
            module_effects.loc[index, "primary_p"]
        )

    score_table = module_scores.copy()
    score_table.insert(0, "matching_ID", design.loc[score_table.index, "matching_ID"])
    score_table.insert(0, "condition", design.loc[score_table.index, "Group"])
    score_table.insert(0, "sample", score_table.index)
    score_long = score_table.melt(
        id_vars=["sample", "condition", "matching_ID"], var_name="module", value_name="score_z"
    )

    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    target_map.to_csv(TABLE_ROOT / "gbs_proteomics_target_aptamer_dictionary.csv", index=False)
    aptamer_effects.to_csv(TABLE_ROOT / "gbs_proteomics_aptamer_effects.csv", index=False)
    gene_scores.rename_axis("sample").reset_index().to_csv(
        TABLE_ROOT / "gbs_proteomics_sample_gene_scores.csv", index=False
    )
    score_long.to_csv(TABLE_ROOT / "gbs_proteomics_sample_module_scores.csv", index=False)
    module_effects.to_csv(TABLE_ROOT / "gbs_proteomics_module_effects.csv", index=False)

    payload = {
        "dataset": "GBS-Proteomics",
        "evidence_tier": "public preprint; longitudinal plasma SomaScan 7K",
        "design_counts": raw["Group"].value_counts().to_dict(),
        "paired_acute_recovery": len(pair_ids),
        "analysis_level": "biological sample; paired by matching_ID for acute versus one-year recovery",
        "module_availability": availability,
        "aptamer_effects": aptamer_effects.to_dict("records"),
        "module_effects": module_effects.to_dict("records"),
        "cautions": [
            "The repository accompanies a preprint and should be interpreted as an external longitudinal validation tier until peer review.",
            "SomaScan measures aptamer binding rather than transcript abundance; multiple aptamers may map to one gene.",
        ],
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUT_JSON}")
    print(raw["Group"].value_counts().to_string())


if __name__ == "__main__":
    main()
