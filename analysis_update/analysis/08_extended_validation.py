#!/usr/bin/env python3
"""Reproducible extensions using only peer-reviewed, publicly reported evidence.

v39 correction: include all 10 paired non-IL8 CSF estimates from Table 4B,
including IL6 and CD5, omitted in v38. IL8 is reported separately in source Table 3.
The Kmezic CSF analysis combines age-, sex- and sample-handling-adjusted
coefficients from independent discovery and replication cohorts. No
participant-level observations are reconstructed. The genetic analysis uses
the deposited donor-level nerve summaries to bootstrap cell localization. The
published nerve evidence records peer-reviewed tissue observations from
PXD056286 and calculates only the binomial confidence interval from 52/55.
"""

from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
OUT = ROOT / "results" / "tables"
RNG = np.random.default_rng(20260904)


def bh(values):
    p = np.asarray(values, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    return out


def csf_meta():
    rows = [
        ("SELE", 1.790, .377, 2.330, .766),
        ("IL2RA", .741, .234, 2.470, .313),
        ("CCL3", .836, .270, .684, .534),
        ("CR2", .692, .273, 1.920, .624),
        ("CD1C", .561, .252, 1.040, .299),
        ("THBD", .656, .317, 1.920, .477),
        ("NRP1", .0965, .0485, .288, .119),
        ("CD38", .635, .352, 1.180, .447),
        ("IL6", .370, .514, 1.890, .427),
        ("CD5", -.0882, .126, 1.510, .514),
    ]
    out = []
    for protein, b1, s1, b2, s2 in rows:
        b = np.array([b1, b2]); se = np.array([s1, s2]); w = 1 / se**2
        est = float(np.sum(w * b) / np.sum(w)); sem = float(np.sqrt(1 / np.sum(w)))
        z = est / sem; p = float(2 * stats.norm.sf(abs(z)))
        qstat = float(np.sum(w * (b - est) ** 2)); df = 1
        i2 = max(0.0, (qstat - df) / qstat * 100) if qstat > 0 else 0.0
        tau2 = max(0.0, (qstat-df)/(np.sum(w)-np.sum(w**2)/np.sum(w)))
        wr = 1/(se**2+tau2); rest = float(np.sum(wr*b)/np.sum(wr)); rse=float(np.sqrt(1/np.sum(wr)))
        out.append(dict(protein=protein, discovery_beta=b1, discovery_se=s1,
                        replication_beta=b2, replication_se=s2,
                        fixed_beta=est, fixed_se=sem, ci_low=est-1.96*sem,
                        ci_high=est+1.96*sem, p_value=p, Q=qstat, I2_percent=i2,
                        tau2_DL=tau2, random_beta_DL=rest,
                        random_ci_low=rest-1.96*rse, random_ci_high=rest+1.96*rse,
                        direction_concordant=(b1>0 and b2>0),
                        inference_level="published adjusted cohort coefficients"))
    d = pd.DataFrame(out)
    d["fdr_BH"] = bh(d.p_value)
    d.to_csv(OUT / "csf_kmezic_adjusted_meta.csv", index=False)
    return d


def genetic_bootstrap(n_boot=5000):
    d = pd.read_csv(SRC / "Genetic_donor_celltype.csv")
    genes = ["CDH4", "DIRAS1", "GNG7", "SLC39A3"]
    d = d[(d.disease == "CIDP") & d.gene.isin(genes)].copy()
    output = []
    for (gene, cell), g in d.groupby(["gene", "cell_group"]):
        vals = g.mean_log2_cp10k_plus1.to_numpy(float)
        det = g.percent_expressing.to_numpy(float)
        idx = RNG.integers(0, len(g), size=(n_boot, len(g)))
        bm = vals[idx].mean(1); bd = det[idx].mean(1)
        output.append(dict(gene=gene, cell_group=cell, n_donors=g["sample"].nunique(),
                           mean_expression=vals.mean(), expression_ci_low=np.quantile(bm,.025),
                           expression_ci_high=np.quantile(bm,.975), mean_percent_expressing=det.mean(),
                           detection_ci_low=np.quantile(bd,.025), detection_ci_high=np.quantile(bd,.975),
                           bootstrap_replicates=n_boot))
    out = pd.DataFrame(output)
    out["expression_rank_within_gene"] = out.groupby("gene").mean_expression.rank(ascending=False, method="min")
    out.to_csv(OUT / "genetic_cell_localization_bootstrap.csv", index=False)
    return out


def external_nerve_validation():
    # Wilson interval, calculated from the published 52/55 observation.
    n, x = 55, 52; z = stats.norm.ppf(.975); phat = x/n
    den = 1 + z*z/n
    ctr = (phat + z*z/(2*n))/den
    half = z*math.sqrt(phat*(1-phat)/n + z*z/(4*n*n))/den
    rows = [
        ("C5b-9 deposition", "Endoneurial capillaries", "52/55 biopsies", phat, ctr-half, ctr+half),
        ("C3 and C6 abundance", "Sural nerve proteome", "Increased in CIDP", np.nan, np.nan, np.nan),
        ("CD68+ macrophages", "Endoneurium", "Cellular infiltration observed", np.nan, np.nan, np.nan),
        ("CD8+ T cells", "Endoneurium", "Cellular infiltration observed", np.nan, np.nan, np.nan),
    ]
    out = pd.DataFrame(rows, columns=["feature","localization","published_result","proportion","ci_low","ci_high"])
    out["source"] = "Stascheit et al. 2025; doi:10.1007/s00401-025-02936-w; PXD056286"
    out["inference_level"] = "published external tissue observations"
    out.to_csv(OUT / "cidp_external_nerve_validation.csv", index=False)
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    c = csf_meta(); g = genetic_bootstrap(); n = external_nerve_validation()
    print(f"CSF meta: {len(c)} proteins; bootstrap: {len(g)} gene-cell rows; published nerve evidence: {len(n)} rows")
