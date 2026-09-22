#!/usr/bin/env python3
"""Exploratory donor level sensitivity analyses. See analysis_plan.md.

Usage: python strengthen_fc.py --source DIR --out DIR
DIR source is the original results directory containing tables and JSON.
No cell is treated as an independent biological replicate.
"""
import argparse, hashlib, itertools, json, platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import scipy

def bh(p):
    p=np.asarray(p,float); q=np.full(p.shape,np.nan); ok=np.isfinite(p)
    ix=np.where(ok)[0]; order=ix[np.argsort(p[ok])]
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return q

def fit(y,d,covariates=()):
    x=[np.ones(len(d)),d.disease.eq('CIDP').to_numpy(float)]
    for c in covariates:
        if c=='center':
            levels=sorted(d.center.unique())
            x.extend(d.center.eq(l).to_numpy(float) for l in levels[1:])
        elif c=='sex': x.append(d.sex.eq('female').to_numpy(float))
        else: x.append(d[c].to_numpy(float))
    X=np.column_stack(x); y=np.asarray(y,float)
    rank=int(np.linalg.matrix_rank(X))
    assert rank==X.shape[1], (covariates,rank,X.shape)
    bread=np.linalg.inv(X.T@X); beta=bread@X.T@y
    h=np.sum((X@bread)*X,axis=1); res=y-X@beta
    vc=bread@(X.T@np.diag((res/(1-h))**2)@X)@bread
    se=float(np.sqrt(max(vc[1,1],0))); df=len(y)-rank
    t=beta[1]/se if se else np.nan
    crit=stats.t.ppf(.975,df)
    return dict(delta=float(beta[1]),se_hc3=se,ci_low=float(beta[1]-crit*se),ci_high=float(beta[1]+crit*se),p_hc3=float(2*stats.t.sf(abs(t),df)),residual_df=df,n=int(len(y)))

def compare(y,d):
    c=d.disease.eq('CIDP').to_numpy(); y=np.asarray(y,float)
    r=fit(y,d);r.update(n_case=int(c.sum()),n_reference=int((~c).sum()),p_mannwhitney=float(stats.mannwhitneyu(y[c],y[~c],alternative='two-sided').pvalue),mean_case=float(y[c].mean()),mean_reference=float(y[~c].mean()))
    return r

def run(source,out):
    out.mkdir(parents=True,exist_ok=True); (out/'private').mkdir(exist_ok=True)
    source=Path(source)
    expr=pd.read_csv(source/'tables/gse285983_sample_targeted_expression.csv')
    mod=pd.read_csv(source/'tables/gse285983_sample_module_scores.csv')
    raw=json.loads((source/'gse285983_targeted_pseudobulk.json').read_text())
    planned_fc=raw['gene_panels']['Fc_receptor']
    fc=[g for g in planned_fc if g in set(expr.gene)]
    e=expr.query('cell_group=="Macrophage"')
    meta=e.drop_duplicates('sample').set_index('sample')[['disease','sex','age','center','n_nuclei']].sort_index()
    mat=e.pivot(index='sample',columns='gene',values='log2_cpm_plus_0_5').loc[meta.index,fc]
    z=(mat-mat.mean())/mat.std(ddof=1).replace(0,np.nan)
    score=z.mean(axis=1)
    previous=mod.query('cell_group=="Macrophage" and panel=="Fc_receptor"').set_index('sample').score_z
    assert np.allclose(score,previous.loc[score.index],atol=1e-10), 'Score reproduction failed'
    take=meta.disease.isin(['CIDP','CIAP']); d=meta.loc[take].copy(); y=score.loc[d.index]
    d['fc_score']=y
    d.reset_index().to_csv(out/'private/fc_donor_analysis.csv',index=False)
    rows=[]
    rows.append(dict(analysis='All donors',**compare(y,d),covariates='none',p_primary=compare(y,d)['p_mannwhitney']))
    dm=d[d.sex.eq('male')];r=compare(dm.fc_score,dm)
    rows.append(dict(analysis='Male donors',**r,covariates='none',p_primary=r['p_mannwhitney']))
    m20=mat.loc[d.index]; z20=(m20-m20.mean())/m20.std(ddof=1).replace(0,np.nan)
    r=compare(z20.mean(axis=1),d)
    rows.append(dict(analysis='Standardisation in CIDP and CIAP only',**r,covariates='none',p_primary=r['p_mannwhitney']))
    for name,cov in [('Age and center',('age','center')),('Age, center and sex',('age','center','sex'))]:
        r=fit(y,d,cov);rows.append(dict(analysis=name,**r,n_case=9,n_reference=11,covariates=', '.join(cov),p_primary=r['p_hc3']))
    groups=[np.flatnonzero(d.center.eq(l)) for l in sorted(d.center.unique())]
    x=d.disease.eq('CIDP').to_numpy(float); yr=np.asarray(y,float).copy(); xr=x.copy()
    for ix in groups: yr[ix]-=yr[ix].mean();xr[ix]-=xr[ix].mean()
    observed=float(xr@yr/(xr@xr)); permutations=[]
    choices=[list(itertools.combinations(ix,int(x[ix].sum()))) for ix in groups]
    for allocations in itertools.product(*choices):
        xp=np.zeros(len(d)); xp[list(itertools.chain.from_iterable(allocations))]=1
        for ix in groups: xp[ix]-=xp[ix].mean()
        permutations.append(float(xp@yr/(xp@xp)))
    pexact=float(np.mean(np.abs(permutations)>=abs(observed)-1e-12))
    rows.append(dict(analysis='Exact permutation within centers',delta=observed,n=20,n_case=9,n_reference=11,covariates='center',p_primary=pexact,permutations=len(permutations)))
    pd.DataFrame(rows).to_csv(out/'fc_sensitivity.csv',index=False)
    genes=[]
    for gene in fc:
        r=compare(mat.loc[d.index,gene],d)
        genes.append(dict(gene=gene,**r,delta_z=float(z.loc[d.index[d.disease.eq('CIDP')],gene].mean()-z.loc[d.index[d.disease.eq('CIAP')],gene].mean()),n_case_above_zero=int((mat.loc[d.index[d.disease.eq('CIDP')],gene]>-1+1e-10).sum()),n_reference_above_zero=int((mat.loc[d.index[d.disease.eq('CIAP')],gene]>-1+1e-10).sum())))
    gt=pd.DataFrame(genes);gt['q_bh_available_genes']=bh(gt.p_mannwhitney);gt.to_csv(out/'fc_component_genes.csv',index=False)
    panels={'Activating receptors':['FCGR1A','FCGR2A','FCGR3A'],'Inhibitory receptor':['FCGR2B'],'Neonatal Fc receptor':['FCGRT'],'Signalling components':['FCER1G','TYROBP','SYK','LYN','HCK']}
    prows=[]
    for panel,gs in panels.items(): prows.append(dict(component=panel,genes=', '.join(gs),**compare(z.loc[d.index,gs].mean(axis=1),d)))
    pt=pd.DataFrame(prows);pt['q_bh_4']=bh(pt.p_mannwhitney);pt.to_csv(out/'fc_functional_components.csv',index=False)
    omitted=[]
    for i,s in enumerate(d.index):
        dd=d.drop(s);omitted.append(dict(type='donor',omitted=f'Anonymous omission {i+1:02d}',**compare(dd.fc_score,dd)))
    for c in sorted(d.center.unique()):
        dd=d[d.center.ne(c)];omitted.append(dict(type='center',omitted=c,**compare(dd.fc_score,dd)))
    for g in fc: omitted.append(dict(type='gene',omitted=g,**compare(z.loc[d.index,[a for a in fc if a!=g]].mean(axis=1),d)))
    omissions=pd.DataFrame(omitted)
    omissions.to_csv(out/'private/fc_omission_analyses.csv',index=False)
    omissions.groupby('type').agg(n_analyses=('delta','size'),delta_min=('delta','min'),delta_max=('delta','max'),p_min=('p_mannwhitney','min'),p_max=('p_mannwhitney','max')).reset_index().to_csv(out/'fc_omission_summary.csv',index=False)
    omissions[omissions.type.ne('donor')].to_csv(out/'fc_gene_center_omissions.csv',index=False)
    broad=pd.DataFrame(raw['module_comparisons']);broad=broad[broad.comparison.eq('CIDP_vs_CIAP')].copy()
    broad['q_bh_all_cell_module_tests']=bh(broad.p_value)
    broad.to_csv(out/'module_global_multiplicity.csv',index=False)
    dem=[]
    for dis,ds in d.groupby('disease'):
        dem.append(dict(disease=dis,characteristic='Donors',value=len(ds)))
        dem.extend(dict(disease=dis,characteristic=f'Sex: {s}',value=int(ds.sex.eq(s).sum())) for s in ['male','female'])
        dem.extend(dict(disease=dis,characteristic=f'Center: {c}',value=int(ds.center.eq(c).sum())) for c in sorted(d.center.unique()))
        dem.extend([dict(disease=dis,characteristic='Age mean in years',value=float(ds.age.mean())),dict(disease=dis,characteristic='Macrophage nuclei minimum',value=int(ds.n_nuclei.min())),dict(disease=dis,characteristic='Macrophage nuclei maximum',value=int(ds.n_nuclei.max()))])
    pd.DataFrame(dem).to_csv(out/'design_balance.csv',index=False)
    summary={'analysis_date':'2026-09-07','exploratory':True,'original_score_max_absolute_error':float(np.max(np.abs(score-previous.loc[score.index]))),'n_genes':len(fc),'n_all_diagnoses':len(meta),'original_family':'8 modules within each cell group and contrast','global_family_n':len(broad),'fc_global_q':float(broad.query('cell_group=="Macrophage" and panel=="Fc_receptor"').q_bh_all_cell_module_tests.iloc[0]),'sensitivity':rows,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source/'tables/gse285983_sample_targeted_expression.csv',source/'tables/gse285983_sample_module_scores.csv',source/'gse285983_targeted_pseudobulk.json']}}
    (out/'fc_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps(summary,indent=2,allow_nan=False))
    print(gt[['gene','delta','p_mannwhitney','q_bh_available_genes']].to_string(index=False))
    print(pt[['component','delta','p_mannwhitney','q_bh_4']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.source,a.out)
