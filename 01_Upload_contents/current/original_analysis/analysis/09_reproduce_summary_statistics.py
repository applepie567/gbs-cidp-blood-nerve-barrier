#!/usr/bin/env python3
"""Recalculate v36 summary statistics from the archived cohort/donor tables.

This is a table-level reproduction. It does not reconstruct unavailable raw
matrices, fit new CSF participant models, or perform an independent GWAS.
"""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
import pandas as pd
from scipy import stats
from common import bh_adjust

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'source_data'
OUT=ROOT/'results/tables'

def meta(g,v):
    k=len(g);w=1/v;m=np.average(g,weights=w);q=np.sum(w*(g-m)**2)
    tau=max(0,(q-k+1)/(w.sum()-(w*w).sum()/w.sum()))
    wr=1/(v+tau);est=np.average(g,weights=wr)
    # Modified Hartung-Knapp: residual scale is bounded below by one.
    se=np.sqrt(max(1,np.sum(wr*(g-est)**2)/(k-1))/wr.sum())
    ci=stats.t.ppf(.975,k-1)*se
    return {'estimate':est,'low':est-ci,'high':est+ci,
            'p':2*stats.t.sf(abs(est/se),k-1),'tau2':tau,
            'I2':max(0,(q-k+1)/q)*100 if q>0 else 0}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    b=pd.read_csv(SRC/'Blood_cohort_effects.csv')
    expected=pd.read_csv(SRC/'Blood_meta.csv').set_index('Program')
    rows=[];loco=[]
    for program,d in b.groupby('Program',sort=False):
        m=meta(d['Hedges g'].to_numpy(),d.Variance.to_numpy())
        for key,col in [('estimate','Summary Hedges g'),('low','95% CI low'),('high','95% CI high'),('p','Hartung–Knapp P'),('I2','I² (%)')]:
            np.testing.assert_allclose(m[key],expected.loc[program,col],atol=1e-6,rtol=1e-6)
        rows.append({'Program':program,**m,'positive_cohorts':int((d['Hedges g']>0).sum())})
        for cohort in d.Cohort:
            sub=d[d.Cohort!=cohort];est=meta(sub['Hedges g'].to_numpy(),sub.Variance.to_numpy())['estimate']
            loco.append({'Program':program,'Omitted cohort':cohort,'Summary Hedges g':est})
    pd.DataFrame(rows).to_csv(OUT/'blood_meta_recomputed.csv',index=False)
    pd.DataFrame(loco).to_csv(SRC/'Figure_2C_leave_one_cohort_out.csv',index=False)

    nerve=pd.read_csv(OUT/'gse285983_sample_module_scores.csv')
    rows=[]
    for (cell,panel),d in nerve.groupby(['cell_group','panel'],sort=False):
        x=d.loc[d.disease=='CIDP','score_z'].dropna();y=d.loc[d.disease=='CIAP','score_z'].dropna()
        if len(x)<2 or len(y)<2:continue
        rows.append(dict(cell_group=cell,panel=panel,n_case=len(x),n_reference=len(y),
            delta_mean_z=x.mean()-y.mean(),p_value=stats.mannwhitneyu(x,y,alternative='two-sided').pvalue))
    n=pd.DataFrame(rows)
    n['fdr_within_celltype_comparison']=n.groupby('cell_group').p_value.transform(lambda p:bh_adjust(p))
    old=pd.read_csv(SRC/'CIDP_module_effects.csv').query("comparison=='CIDP_vs_CIAP'")
    joined=n.merge(old,on=['cell_group','panel'],suffixes=('_new','_old'))
    for col in ['delta_mean_z','p_value','fdr_within_celltype_comparison']:
        np.testing.assert_allclose(joined[col+'_new'],joined[col+'_old'],atol=1e-10)
    n.to_csv(OUT/'cidp_module_effects_recomputed.csv',index=False)

    donor=pd.read_csv(SRC/'CIDP_donor_heterogeneity.csv').set_index('sample')
    rows=[]
    for a,b in combinations(donor.columns,2):
        d=donor[[a,b]].dropna();r,p=stats.spearmanr(d[a],d[b])
        rows.append(dict(module_1=a,module_2=b,n_donors=len(d),rho=r,p_value=p))
    c=pd.DataFrame(rows);c['fdr_BH_10_pairs']=bh_adjust(c.p_value)
    c.to_csv(OUT/'cidp_donor_module_correlations.csv',index=False)
    cl=pd.read_csv(SRC/'CIDP_clinical.csv');cl['fdr_BH_15_tests']=bh_adjust(cl.p_value)
    cl.to_csv(OUT/'cidp_clinical_correlations_with_fdr.csv',index=False)
    coverage=pd.read_csv(SRC/'CSF_published_evidence.csv').groupby('program').agg(
        studies=('study','nunique'),platforms=('platform','nunique')).loc[
            ['Protein accumulation','Structural proteolysis','Local inflammation',
             'Complement','Protein handling and neural structure']]
    coverage.to_csv(SRC/'Figure_3C_coverage.csv')
    report={'level':'archived cohort and donor summaries','raw_matrix_pipeline_rerun':False,
        'blood_meta_modules_matched':len(expected),
        'nerve_contrasts_matched':len(joined),'donor_pairs':len(c),'clinical_tests':len(cl),
        'donor_pair_FDR_significant':int((c.fdr_BH_10_pairs<.05).sum()),
        'clinical_FDR_significant':int((cl.fdr_BH_15_tests<.05).sum())}
    (ROOT/'metadata/SUMMARY_REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
