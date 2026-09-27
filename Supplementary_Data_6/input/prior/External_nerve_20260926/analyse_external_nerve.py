#!/usr/bin/env python3
"""Frozen 11-gene Fc analysis in GSE213455. See analysis_plan.md.

Run: python analyse_external_nerve.py --root .
Dependencies: numpy, pandas, scipy, matplotlib (no statsmodels required).
"""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats, special

GENES = ['FCGR1A', 'FCGR2A', 'FCGR2B', 'FCGR3A', 'FCGR3B', 'FCGRT',
         'FCER1G', 'TYROBP', 'SYK', 'LYN', 'HCK']


def hc3(y, disease, covariates=None):
    """OLS disease coefficient and HC3 sandwich CI, using residual t df."""
    X = np.column_stack([np.ones(len(y)), disease] + ([] if covariates is None else list(covariates)))
    assert np.linalg.matrix_rank(X) == X.shape[1]
    bread = np.linalg.inv(X.T @ X)
    beta = bread @ X.T @ y
    resid = y - X @ beta
    leverage = np.einsum('ij,jk,ik->i', X, bread, X)
    cov = bread @ (X.T @ (X * (resid / (1 - leverage))[:, None] ** 2)) @ bread
    se = np.sqrt(cov[1, 1])
    df = len(y) - X.shape[1]
    critical = stats.t.ppf(0.975, df)
    return dict(difference=float(beta[1]), hc3_se=float(se),
                ci_low=float(beta[1] - critical * se), ci_high=float(beta[1] + critical * se),
                model_p=float(2 * stats.t.sf(abs(beta[1] / se), df)), residual_df=int(df))


def exact_rank_p(y, disease):
    """Exhaustive two-sided Mann–Whitney rank permutation, including ties."""
    ranks = stats.rankdata(y)
    n1, n = int(disease.sum()), len(y)
    expected = n1 * (n + 1) / 2
    observed = abs(ranks[disease == 1].sum() - expected)
    possible = np.array([abs(ranks[list(c)].sum() - expected)
                         for c in itertools.combinations(range(n), n1)])
    return float((possible >= observed - 1e-12).mean()), len(possible)


def hedges_g(y, disease):
    a, b = y[disease == 1], y[disease == 0]
    df = len(a) + len(b) - 2
    pooled = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / df)
    correction = np.exp(special.gammaln(df/2) - .5*np.log(df/2) - special.gammaln((df-1)/2))
    return float(correction * (a.mean() - b.mean()) / pooled)


def bh(values):
    p = np.array(values)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted,1)
    return out


def score_frame(frame):
    return ((frame - frame.mean()) / frame.std(ddof=1)).mean(axis=1)


def parse_samples(path):
    lines=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'): break
            if line.startswith('!Sample_'): lines.append(next(csv.reader([line],delimiter='\t')))
    fields={}
    for row in lines:
        if row[0]=='!Sample_characteristics_ch1':
            key = row[1].split(': ',1)[0]
            fields[key]=[v.split(': ',1)[1] for v in row[1:]]
        elif row[0] in ['!Sample_geo_accession','!Sample_title']:
            fields[row[0].removeprefix('!Sample_')]=row[1:]
    meta=pd.DataFrame(fields).rename(columns={'geo_accession':'sample_id','pathological diagnosis':'deposited_diagnosis','gender':'sex'})
    mapping={'vasculitic neuropathy':'VN','inflammatory neuropathy':'CIDP','Lymphomatous neuropathy':'NL'}
    meta['diagnosis']=meta.deposited_diagnosis.map(mapping)
    assert meta.diagnosis.notna().all()
    meta['mapping_evidence']='GEO sample labels and Cerri et al. 2022 Methods (doi:10.3389/fonc.2022.974751)'
    meta['primary_eligible']=meta.diagnosis.isin(['CIDP','VN'])
    meta['exclusion_reason']=np.where(meta.primary_eligible,'','NL excluded by frozen analysis plan')
    meta['age']=meta.age.astype(int)
    return meta


def main(root):
    out=root/'results'
    out.mkdir(parents=True,exist_ok=True)
    matrix_path=root/'inputs/GSE213455_series_matrix.txt.gz'
    ann_path=root/'inputs/GPL13369_annotation.tsv'
    if not ann_path.exists(): ann_path=ann_path.with_suffix('.tsv.gz')
    meta=parse_samples(matrix_path)
    assert meta.groupby('diagnosis').size().to_dict()=={'CIDP':4,'NL':3,'VN':9}
    meta.to_csv(out/'GSE213455_sample_mapping.csv',index=False)
    eligible=meta[meta.primary_eligible].set_index('sample_id')
    matrix=pd.read_csv(matrix_path,sep='\t',comment='!',index_col=0)
    matrix.index=matrix.index.astype(str)
    assert list(matrix.columns)==list(meta.sample_id)
    assert not matrix.index.duplicated().any()
    ann=pd.read_csv(ann_path,sep='\t',dtype=str,low_memory=False).fillna('')
    # Only exact single-symbol matches are used. No effect-based probe selection.
    probes=ann[ann.SYMBOL.isin(GENES)].copy()
    probes['present_in_matrix']=probes.ID.isin(matrix.index)
    probes['mapping_rule']='exact single gene symbol from GPL13369'
    probes.to_csv(out/'Fc_probe_annotation.csv',index=False)
    probes=probes[probes.present_in_matrix]
    assert not probes.ID.duplicated().any()
    median,mean,coverage={},{},[]
    for gene in GENES:
        ids=probes.loc[probes.SYMBOL==gene,'ID'].tolist()
        values=matrix.loc[ids,eligible.index]
        valid=len(ids)>0 and np.isfinite(values.to_numpy()).all()
        gene_sd=np.nan
        if valid:
            med=values.median(axis=0)
            gene_sd=med.std(ddof=1)
            valid=bool(gene_sd>0)
            if valid:
                median[gene]=med
                mean[gene]=values.mean(axis=0)
        coverage.append(dict(gene=gene,n_probes=len(ids),probe_ids=';'.join(ids),included=valid,
                             sd_median_expression=gene_sd,
                             reason='included' if valid else 'absent, missing or invariant'))
    coverage=pd.DataFrame(coverage)
    coverage.to_csv(out/'Fc_gene_coverage.csv',index=False)
    assert len(median)>=8, 'Coverage below threshold; module test not permitted'
    expr=pd.DataFrame(median)
    expr.index.name='sample_id'
    expr.to_csv(out/'Fc_gene_expression_median.csv')
    pd.DataFrame(mean).to_csv(out/'Fc_gene_expression_mean.csv')
    z=(expr-expr.mean())/expr.std(ddof=1)
    z.to_csv(out/'Fc_gene_zscores.csv')
    disease=(eligible.diagnosis=='CIDP').to_numpy().astype(int)
    scores=score_frame(expr).to_numpy()
    alternative=score_frame(pd.DataFrame(mean)).to_numpy()
    donors=eligible[['diagnosis','age','sex']].copy()
    donors['Fc_score']=scores
    donors['Fc_score_mean_probe']=alternative
    donors.to_csv(out/'Fc_donor_scores.csv')
    primary=hc3(scores,disease)
    p,nperm=exact_rank_p(scores,disease)
    primary.update(model='Unadjusted (median probes)',p_exact_rank=p,permutations=nperm,
                   hedges_g=hedges_g(scores,disease),n_CIDP=4,n_VN=9,n_genes=len(median))
    adjusted=hc3(scores,disease,[(eligible.age-eligible.age.mean()).to_numpy(),(eligible.sex=='M').to_numpy().astype(int)])
    adjusted.update(model='Age and sex adjusted (median probes)',p_exact_rank=np.nan,hedges_g=np.nan)
    alt=hc3(alternative,disease)
    alt.update(model='Unadjusted (mean probes)',p_exact_rank=exact_rank_p(alternative,disease)[0],hedges_g=hedges_g(alternative,disease))
    pd.DataFrame([primary,adjusted,alt]).to_csv(out/'Fc_module_results.csv',index=False)
    genes=[]
    for gene in expr:
        y=expr[gene].to_numpy()
        r=hc3(y,disease)
        r.update(gene=gene,n_CIDP=4,n_VN=9,hedges_g=hedges_g(y,disease),p_exact_rank=exact_rank_p(y,disease)[0])
        genes.append(r)
    gene_results=pd.DataFrame(genes)
    gene_results['q_BH_11']=bh(gene_results.p_exact_rank)
    gene_results.to_csv(out/'Fc_gene_results.csv',index=False)
    # Recompute standardisation after each omission, without selecting genes by outcome.
    omissions=[]
    for i,sample in enumerate(expr.index):
        retained=np.arange(len(expr))!=i
        s=score_frame(expr.iloc[retained]).to_numpy()
        d=disease[retained]
        omissions.append(dict(omitted_sample=sample,omitted_diagnosis=eligible.iloc[i].diagnosis,
                              difference=float(s[d==1].mean()-s[d==0].mean()),hedges_g=hedges_g(s,d)))
    loo=pd.DataFrame(omissions)
    loo.to_csv(out/'Fc_leave_one_donor_out.csv',index=False)
    qc=dict(n_matrix_probes=len(matrix),n_matrix_samples=matrix.shape[1],missing_values=int(matrix.isna().sum().sum()),
            eligible_samples=len(eligible),covered_genes=list(expr.columns),n_unique_sample_ids=meta.sample_id.nunique(),
            value_min=float(matrix.min().min()),value_max=float(matrix.max().max()),
            unique_quantile_profiles=int(matrix.apply(lambda x: tuple(np.sort(x))).T.drop_duplicates().shape[0]),
            score_z_gene_means_abs_max=float(z.mean().abs().max()),
            input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [matrix_path,ann_path]},
            analysis_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            primary=primary,adjusted=adjusted,alternative=alt,
            loo_min=float(loo.difference.min()),loo_max=float(loo.difference.max()),
            direction='CIDP minus vasculitic neuropathy',interpretation_scope='independent whole-nerve comparator extension; not cell-specific CIDP-versus-CIAP replication')
    # JSON null denotes not applicable, not a computed zero.
    clean=json.loads(json.dumps(qc),parse_constant=lambda _:None)
    (out/'analysis_summary.json').write_text(json.dumps(clean,indent=2,allow_nan=False),encoding='utf8')
    print(json.dumps(clean,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    main(parser.parse_args().root)
