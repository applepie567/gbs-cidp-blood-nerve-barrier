#!/usr/bin/env python3
"""BNB identity anchor, rat EAN time course, and published-table validation.

GSE107574 is a normal human endoneurial endothelial/microvessel transcriptome
and is used only to anchor blood-nerve-barrier identity. It is not a disease
contrast and is never labelled as BBB. GSE133750 is a rat sciatic-nerve EAN
time course and is retained as an animal mechanistic tier. OEP002315/OEP002701
is represented by author-supplied supplementary GSEA tables, not by a claimed
raw-data reanalysis.
"""

from __future__ import annotations

import gzip
import json
import shutil
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from common import GENE_MODULES, bh_adjust, unpaired_effect  # noqa: E402


TABLE_ROOT = ROOT / "results/tables"
OUT_JSON = ROOT / "results/bnb_ean_external_validation.json"


def read_gzipped_xlsx(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rb") as source, tempfile.NamedTemporaryFile(suffix=".xlsx") as target:
        shutil.copyfileobj(source, target)
        target.flush()
        return pd.read_excel(target.name)


def analyze_bnb_anchor() -> dict:
    path = ROOT / "data/raw/GSE107574/GSE107574_fpkm_gene_expression_values.xlsx.gz"
    frame = read_gzipped_xlsx(path)
    frame["gene"] = frame["gene_id"].astype(str).str.rsplit("_", n=1).str[-1]
    sample_groups = {
        "P3pHEnd_EC": "cultured_endoneurial_endothelial_cell",
        "P8_pHEnd_EC_basal": "cultured_endoneurial_endothelial_cell",
        "32P1": "LCM_endoneurial_microvessel",
        "203P1": "LCM_endoneurial_microvessel",
        "346P1": "LCM_endoneurial_microvessel",
        "347P1": "LCM_endoneurial_microvessel",
    }
    target_genes = sorted({gene for genes in GENE_MODULES.values() for gene in genes})
    expression = frame.loc[frame["gene"].isin(target_genes), ["gene", *sample_groups]].copy()
    expression = expression.groupby("gene", sort=False).sum(numeric_only=True)
    long = expression.rename_axis("gene").reset_index().melt(
        id_vars="gene", var_name="sample", value_name="fpkm"
    )
    long["preparation"] = long["sample"].map(sample_groups)
    long["modules"] = long["gene"].map(
        lambda gene: ";".join(module for module, genes in GENE_MODULES.items() if gene in genes)
    )
    summary = long.groupby(["gene", "preparation"], as_index=False).agg(
        median_fpkm=("fpkm", "median"),
        mean_fpkm=("fpkm", "mean"),
        detection_fraction=("fpkm", lambda values: float(np.mean(np.asarray(values) > 0))),
        n_preparations=("fpkm", "size"),
    )
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    long.to_csv(TABLE_ROOT / "gse107574_bnb_target_expression.csv", index=False)
    summary.to_csv(TABLE_ROOT / "gse107574_bnb_target_summary.csv", index=False)
    return {
        "dataset": "GSE107574",
        "role": "normal human BNB endothelial identity anchor; descriptive only",
        "design_counts": {"cultured_endoneurial_endothelial_cell": 2, "LCM_endoneurial_microvessel": 4},
        "no_disease_contrast": True,
        "not_bbb": True,
        "target_summary": summary.to_dict("records"),
    }


def query_rat_map(symbols: list[str]) -> pd.DataFrame:
    cache = ROOT / "data/raw/GSE133750/rat_target_gene_mapping.csv"
    if cache.exists():
        return pd.read_csv(cache)
    rows: list[dict] = []
    for start in range(0, len(symbols), 30):
        chunk = symbols[start:start + 30]
        query = "symbol:(" + " OR ".join(chunk) + ")"
        params = urllib.parse.urlencode({
            "q": query,
            "species": "10116",
            "fields": "symbol,ensembl.gene",
            "size": "100",
        })
        with urllib.request.urlopen("https://mygene.info/v3/query?" + params, timeout=60) as response:
            payload = json.load(response)
        for hit in payload.get("hits", []):
            ensembl = hit.get("ensembl")
            if isinstance(ensembl, list):
                ensembl = ensembl[0] if ensembl else {}
            gene_id = ensembl.get("gene", "") if isinstance(ensembl, dict) else ""
            if gene_id:
                rows.append({"rat_symbol": hit.get("symbol", ""), "rat_ensembl": gene_id})
    found = pd.DataFrame(rows).drop_duplicates()
    mapped = []
    for symbol in symbols:
        exact = found.loc[found["rat_symbol"].str.upper().eq(symbol.upper())]
        if not exact.empty:
            for row in exact.itertuples(index=False):
                mapped.append({"query_symbol": symbol, "rat_symbol": row.rat_symbol, "rat_ensembl": row.rat_ensembl})
    out = pd.DataFrame(mapped).drop_duplicates()
    out.to_csv(cache, index=False)
    return out


def analyze_ean() -> dict:
    path = ROOT / "data/raw/GSE133750/GSE133750_fpkm.txt.gz"
    fpkm = pd.read_csv(path, sep="\t", index_col=0)
    modules = dict(GENE_MODULES)
    # Rats do not have a direct CXCL8 orthologue; ELR+ chemokines Cxcl1/Cxcl2
    # provide the prespecified functional surrogate for this animal-only tier.
    modules["CXCL8_CXCR1_2"] = ["Cxcl1", "Cxcl2", "Cxcr1", "Cxcr2"]
    query_symbols = sorted({gene for genes in modules.values() for gene in genes})
    mapping = query_rat_map(query_symbols)
    mapping = mapping.loc[mapping["rat_ensembl"].isin(fpkm.index)].copy()
    mapping.to_csv(TABLE_ROOT / "gse133750_rat_target_gene_mapping_used.csv", index=False)
    log_expression = np.log2(fpkm + 0.5)
    z = log_expression.sub(log_expression.mean(axis=1), axis=0).div(log_expression.std(axis=1, ddof=1), axis=0)
    groups = {
        "Control": [column for column in fpkm.columns if column.startswith("C")],
        "Early_neuritis": [column for column in fpkm.columns if column.startswith("E")],
        "Peak_neuritis": [column for column in fpkm.columns if column.startswith("P")],
        "Late_neuritis": [column for column in fpkm.columns if column.startswith("L")],
    }
    comparisons = [
        ("Early_neuritis_vs_Control", "Early_neuritis"),
        ("Peak_neuritis_vs_Control", "Peak_neuritis"),
        ("Late_neuritis_vs_Control", "Late_neuritis"),
    ]

    gene_rows: list[dict] = []
    for map_row in mapping.itertuples(index=False):
        values = log_expression.loc[map_row.rat_ensembl]
        for comparison, case_group in comparisons:
            gene_rows.append({
                "dataset": "GSE133750",
                "comparison": comparison,
                "query_symbol": map_row.query_symbol,
                "rat_symbol": map_row.rat_symbol,
                "rat_ensembl": map_row.rat_ensembl,
                "modules": ";".join(module for module, genes in modules.items() if map_row.query_symbol in genes),
                **unpaired_effect(values.loc[groups[case_group]], values.loc[groups["Control"]]),
            })
    gene_effects = pd.DataFrame(gene_rows)
    for comparison, index in gene_effects.groupby("comparison").groups.items():
        gene_effects.loc[index, "targeted_fdr_within_comparison"] = bh_adjust(
            gene_effects.loc[index, "permutation_p"]
        )

    module_scores = {}
    availability = {}
    for module, symbols in modules.items():
        ids = mapping.loc[mapping["query_symbol"].isin(symbols), "rat_ensembl"].drop_duplicates().tolist()
        availability[module] = mapping.loc[mapping["rat_ensembl"].isin(ids), "rat_symbol"].drop_duplicates().tolist()
        if len(ids) >= 2:
            module_scores[module] = z.loc[ids].mean(axis=0)
    module_scores = pd.DataFrame(module_scores)
    module_rows: list[dict] = []
    for module in module_scores.columns:
        for comparison, case_group in comparisons:
            module_rows.append({
                "dataset": "GSE133750",
                "comparison": comparison,
                "module": module,
                "n_rat_genes": len(availability[module]),
                "rat_genes_available": ";".join(availability[module]),
                **unpaired_effect(
                    module_scores.loc[groups[case_group], module],
                    module_scores.loc[groups["Control"], module],
                ),
            })
    module_effects = pd.DataFrame(module_rows)
    for comparison, index in module_effects.groupby("comparison").groups.items():
        module_effects.loc[index, "module_fdr_within_comparison"] = bh_adjust(
            module_effects.loc[index, "permutation_p"]
        )
    score_long = module_scores.rename_axis("sample").reset_index().melt(
        id_vars="sample", var_name="module", value_name="score_z"
    )
    score_long["condition"] = score_long["sample"].map(
        {sample: group for group, samples in groups.items() for sample in samples}
    )
    gene_effects.to_csv(TABLE_ROOT / "gse133750_gene_effects.csv", index=False)
    module_effects.to_csv(TABLE_ROOT / "gse133750_module_effects.csv", index=False)
    score_long.to_csv(TABLE_ROOT / "gse133750_sample_module_scores.csv", index=False)
    return {
        "dataset": "GSE133750",
        "role": "rat EAN sciatic-nerve mechanistic time course",
        "design_counts": {group: len(samples) for group, samples in groups.items()},
        "cxcl8_surrogate": "Cxcl1/Cxcl2-Cxcr1/Cxcr2; rat has no direct CXCL8 orthologue",
        "module_availability": availability,
        "gene_effects": gene_effects.to_dict("records"),
        "module_effects": module_effects.to_dict("records"),
    }


def extract_oep_tables() -> dict:
    base = ROOT / "data/raw/OEP002315_OEP002701/supplement"
    gsea_path = base / "41598_2023_32427_MOESM6_ESM.xls"
    selected_rows = []
    keywords = r"complement|coagulation|fc gamma|leukocyte transendothelial|chemokine|phagosome|endocytosis"
    workbook = pd.ExcelFile(gsea_path)
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(gsea_path, sheet_name=sheet)
        hits = frame.loc[frame["Description"].astype(str).str.contains(keywords, case=False, regex=True)].copy()
        for _, row in hits.iterrows():
            selected_rows.append({
                "dataset": "OEP002315/OEP002701",
                "cell_subset": sheet,
                "pathway_id": row["ID"],
                "description": row["Description"],
                "NES": float(row["NES"]),
                "pvalue": float(row["pvalue"]),
                "p_adjust": float(row["p.adjust"]),
                "source_level": "author-supplied supplementary GSEA result",
            })
    selected = pd.DataFrame(selected_rows).sort_values(["cell_subset", "p_adjust", "description"])
    selected.to_csv(TABLE_ROOT / "oep002315_oep002701_selected_author_gsea.csv", index=False)

    clinical_path = base / "41598_2023_32427_MOESM4_ESM.xlsx"
    clinical = pd.read_excel(clinical_path, sheet_name="Clinical Characteristics")
    clinical = clinical.rename(columns={clinical.columns[0]: "donor"})
    clinical = clinical.loc[clinical["donor"].astype(str).str.match(r"^(T|HC)\d+$")].copy()
    clinical.to_csv(TABLE_ROOT / "oep002315_oep002701_clinical_metadata.csv", index=False)
    return {
        "dataset": "OEP002315/OEP002701",
        "role": "published monocyte GSEA external validation; no raw matrix reanalysis claimed",
        "design": "3 peak AIDP, 2 late/convalescent AIDP, 3 healthy controls",
        "selected_author_gsea": selected.to_dict("records"),
    }


def write_evidence_hierarchy() -> pd.DataFrame:
    rows = [
        ["GSE304871", "Human GBS", "sorted CD11b/CD4/CD8 bulk RNA-seq", "early untreated AIDP vs HC", "Direct reanalysis", "Same-study counterpart GSE304872 is not independent"],
        ["PRJNA1293757", "Human GBS", "PBMC scRNA-seq", "3 treatment-naive acute AIDP vs 2 HC", "Direct sample-level pseudobulk", "Exploratory; exact P values are coarse"],
        ["GSE211225", "Human GBS", "whole-blood bulk RNA-seq", "6 acute, 10 post-acute, 6 HC", "Direct reanalysis", "Cross-sectional; acute/post-acute are not paired"],
        ["GSE31014", "Human GBS", "peripheral-leukocyte microarray", "7 GBS vs 7 HC", "Direct reanalysis", "Processed centered matrix; cross-platform directional validation"],
        ["GBS-Proteomics", "Human GBS", "plasma SomaScan 7K", "20 acute + same 20 at 1 year + 15 HC", "Direct paired targeted reanalysis", "Public preprint tier; aptamer abundance"],
        ["GSE285983", "Human CIDP/PNS", "sural-nerve snRNA-seq", "37 donors including 9 CIDP", "Direct donor-level pseudobulk", "Tissue/diagnosis heterogeneity; not longitudinal IVIg"],
        ["GSE107574", "Normal human BNB", "endoneurial EC and LCM microvessels", "2 cultured EC + 4 microvessels", "Direct descriptive identity anchor", "Not a disease contrast and not BBB"],
        ["OEP002315/OEP002701", "Human AIDP", "PBMC scRNA-seq supplementary GSEA", "3 peak, 2 late, 3 HC", "Author-result-table validation", "No individual expression matrix included in this package"],
        ["KSX10000018-21", "Human CIDP", "2026 stable-treated CIDP multi-omic study", "20 CIDP vs 20 HC; scRNA subset 2+2", "Publication-level external validation", "Stable on maintenance immunotherapy; not active/untreated/longitudinal IVIg"],
        ["GSE133750", "Rat EAN", "sciatic-nerve bulk RNA-seq time course", "control/early/peak/late, n=3 each", "Direct animal mechanistic reanalysis", "Animal sensitivity tier; CXCL8 represented by ELR+ surrogate"],
        ["GSE304872", "Human GBS", "PBMC scRNA-seq Seurat object", "same donors/study as GSE304871", "Same-study published/derived support", "Not counted as independent validation"],
        ["CIDP managed-access/IVIg cohorts", "Human CIDP", "controlled/managed individual-level data", "access pending", "Access gate", "No treatment-response claim until individual data are obtained"],
    ]
    frame = pd.DataFrame(rows, columns=["resource", "biological_scope", "modality", "design", "evidence_tier", "guardrail"])
    frame.to_csv(TABLE_ROOT / "upgraded_evidence_hierarchy.csv", index=False)
    return frame


def main() -> None:
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        "bnb_anchor": analyze_bnb_anchor(),
        "ean": analyze_ean(),
        "oep_external": extract_oep_tables(),
        "evidence_hierarchy": write_evidence_hierarchy().to_dict("records"),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
