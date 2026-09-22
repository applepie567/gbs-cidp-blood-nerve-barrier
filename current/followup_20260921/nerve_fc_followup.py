"""Focused reference comparison and direct standardisation feasibility audit."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import numpy as np
import pandas as pd
from scipy import stats
import scipy

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'archive_original/GBS_CIDP_v41_author_archive_2026-09-07/03_private_reanalysis_workspace'
EXT=ROOT/'extension_original/GBS_CIDP_exploratory_extension_20260921'
OUT=ROOT/'results'
PRIVATE=ROOT/'private'
FC=['FCGR1A','FCGR2A','FCGR2B','FCGR3A','FCGR3B','FCGRT','FCER1G','TYROBP','SYK','LYN','HCK']


def bh(p):
    p=np.asarray(p,float)
    ix=np.argsort(p)
    q=np.empty(len(p))
    q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q


def bootstrap_difference(a,b,rng,n=20000):
    return np.quantile(a[rng.integers(len(a),size=(n,len(a)))].mean(axis=1)-
                       b[rng.integers(len(b),size=(n,len(b)))].mean(axis=1),[.025,.975])


def run():
    OUT.mkdir(exist_ok=True);PRIVATE.mkdir(exist_ok=True)
    ex=pd.read_csv(SOURCE/'original_results/tables/gse285983_sample_targeted_expression.csv')
    ex=ex[ex.cell_group.eq('Macrophage') & ex.gene.isin(FC)]
    assert not ex.duplicated(['sample','gene']).any()
    assert ex.n_nuclei.min()>=20
    m=ex.pivot(index='sample',columns='gene',values='log2_cpm_plus_0_5')[FC]
    assert m.shape==(37,11) and m.notna().all().all()
    d=ex.drop_duplicates('sample').set_index('sample').loc[m.index]
    z=(m-m.mean())/m.std(ddof=1)
    d['score']=z.mean(axis=1)
    archived=pd.read_csv(SOURCE/'original_results/tables/gse285983_sample_module_scores.csv')
    archived=archived[archived.cell_group.eq('Macrophage') & archived.panel.eq('Fc_receptor')].set_index('sample')
    error=float(np.max(np.abs(d.score-archived.loc[d.index,'score_z'])))
    assert error<1e-12
    focus=d[d.disease.isin(['CIDP','CIAP','CTRL'])].copy()
    assert focus.disease.value_counts().to_dict()=={'CIAP':11,'CIDP':9,'CTRL':4}
    focus[['disease','sex','age','center','n_nuclei','score']].to_csv(PRIVATE/'fc_reference_donor_scores.csv')
    rng=np.random.default_rng(20260921)
    rows=[]
    for case,ref in [('CIDP','CIAP'),('CIDP','CTRL'),('CIAP','CTRL')]:
        a=focus.loc[focus.disease.eq(case),'score'].to_numpy()
        b=focus.loc[focus.disease.eq(ref),'score'].to_numpy()
        assert len(np.unique(np.r_[a,b]))==len(a)+len(b), 'Ties need permutation implementation'
        ci=bootstrap_difference(a,b,rng)
        rows.append(dict(case=case,reference=ref,n_case=len(a),n_reference=len(b),
                         mean_case=a.mean(),mean_reference=b.mean(),delta=a.mean()-b.mean(),
                         bootstrap_low=ci[0],bootstrap_high=ci[1],
                         p_mannwhitney_exact=stats.mannwhitneyu(a,b,alternative='two-sided',method='exact').pvalue,
                         p_mannwhitney_original_method=stats.mannwhitneyu(a,b,alternative='two-sided').pvalue,
                         role='Prior contrast' if ref=='CIAP' else 'Exploratory graft reference comparison'))
    result=pd.DataFrame(rows)
    prior=pd.read_csv(EXT/'legacy_reference/strengthening_v41/results/fc_sensitivity.csv')
    baseline=prior.iloc[0]
    assert abs(result.iloc[0].delta+0.673)<.001
    assert abs(result.iloc[0].p_mannwhitney_original_method-0.001087436576734)<1e-12
    result['q_bh_two_reference_tests']=np.nan
    result.loc[1:,'q_bh_two_reference_tests']=bh(result.loc[1:,'p_mannwhitney_exact'])
    old=pd.read_csv(EXT/'results/original_nerve_multiplicity_89.csv')
    new=pd.read_csv(EXT/'results/nerve_program_contrasts.csv')
    assert len(old)==80 and len(new)==9
    p=np.r_[old.p_value,new.p_mannwhitney,result.loc[1:,'p_mannwhitney_exact']]
    assert len(p)==91
    q=bh(p)
    result['q_bh_91']=np.nan
    original_fc=old.index[old.cell_group.eq('Macrophage') & old.panel.eq('Fc_receptor')]
    assert len(original_fc)==1
    result.loc[0,'q_bh_91']=q[original_fc[0]]
    result.loc[1:,'q_bh_91']=q[-2:]
    result.to_csv(OUT/'fc_reference_comparisons.csv',index=False)
    labels=old[['comparison','cell_group','panel']].copy()
    labels['source_family']='Original 80 CIDP versus CIAP module tests'
    extra=new[['cell_group','program']].rename(columns={'program':'panel'}).copy()
    extra['comparison']='CIDP_vs_CIAP'
    extra['source_family']='Nine functional programme extensions'
    refs=pd.DataFrame({'comparison':['CIDP_vs_CTRL','CIAP_vs_CTRL'],
                       'cell_group':['Macrophage']*2,'panel':['Fc_receptor']*2,
                       'source_family':['Focused graft reference comparisons']*2})
    family=pd.concat([labels,extra,refs],ignore_index=True)
    family.insert(0,'test_index',np.arange(1,92))
    family['p_value']=p
    family['q_bh_91']=q
    family.to_csv(OUT/'nerve_multiplicity_91.csv',index=False)
    groups=focus.groupby('disease').agg(n=('score','size'),mean_score=('score','mean'),sd_score=('score','std'),
                                      median_score=('score','median'),min_score=('score','min'),max_score=('score','max'),
                                      mean_age=('age','mean'),n_male=('sex',lambda s:int(s.eq('male').sum())),
                                      nuclei_min=('n_nuclei','min'),nuclei_max=('n_nuclei','max'))
    groups.to_csv(OUT/'fc_reference_group_summary.csv')
    focus.groupby(['disease','center']).size().rename('n_donors').to_csv(OUT/'fc_reference_center_counts.csv')
    cov=pd.read_csv(SOURCE/'results/private/fc_annotation_coverage.csv')
    fine=cov[cov.annotation.str.fullmatch(r'Refined Macro\d+')].copy()
    assert len(fine)==20*18 and not fine.duplicated(['sample','annotation']).any()
    fine['eligible']=fine.n_nuclei.ge(20)
    audit=fine.groupby(['annotation','disease']).agg(n_donors=('sample','size'),eligible_donors=('eligible','sum'),
        min_nuclei=('n_nuclei','min'),max_nuclei=('n_nuclei','max'),zero_donors=('n_nuclei',lambda a:int(a.eq(0).sum())))
    audit.to_csv(OUT/'fc_standardisation_support.csv')
    all_support=fine.groupby('annotation').eligible.sum()
    common=all_support[all_support.eq(20)].index.tolist()
    per_donor=fine.groupby(['sample','disease']).eligible.sum().rename('eligible_subtypes').reset_index()
    per_donor.to_csv(PRIVATE/'fc_donor_subtype_support.csv',index=False)
    support_summary=per_donor.groupby('disease').agg(n_donors=('sample','size'),
        donors_without_eligible_subtype=('eligible_subtypes',lambda v:int(v.eq(0).sum())),
        min_eligible_subtypes=('eligible_subtypes','min'),max_eligible_subtypes=('eligible_subtypes','max'))
    support_summary.to_csv(OUT/'fc_standardisation_group_support.csv')
    summary={'score_reconstruction_max_error':error,'original_donor_count':37,'fc_genes':FC,
             'reference_comparisons':result.replace({np.nan:None}).to_dict('records'),
             'composition_subtypes':18,'common_eligible_subtypes':common,
             'full_cohort_direct_standardisation_estimable':bool(common),
             'standardisation_not_estimated_reason':'No refined macrophage subtype has at least 20 nuclei in every one of the 20 disease donors',
             'seed':20260921,'bootstrap_replicates':20000,'python':platform.python_version(),
             'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,
             'testing_scope':'The 91 tests extend the manuscript family. This is not an omnibus correction for every analysis in the original pipeline.',
             'provenance':'CIDP versus CTRL was calculated in the original pipeline. This follow up reconstructs it and adds CIAP versus CTRL.',
             'input_files':[{'file':str(f.relative_to(SOURCE if f.is_relative_to(SOURCE) else EXT)),
                             'source':'author_archive' if f.is_relative_to(SOURCE) else 'extension',
                             'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
                            for f in [SOURCE/'original_results/tables/gse285983_sample_targeted_expression.csv',
                                      SOURCE/'original_results/tables/gse285983_sample_module_scores.csv',
                                      SOURCE/'results/private/fc_annotation_coverage.csv',
                                      EXT/'legacy_reference/strengthening_v41/results/fc_sensitivity.csv',
                                      EXT/'results/original_nerve_multiplicity_89.csv',
                                      EXT/'results/nerve_program_contrasts.csv']]}
    (OUT/'nerve_followup_summary.json').write_text(json.dumps(summary,indent=2))
    print(result.to_string(index=False))
    print(groups.to_string())
    print('Common subtypes eligible in all 20 donors:',common)
    print(support_summary.to_string())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=SOURCE,help='03_private_reanalysis_workspace in the v41 author archive')
    ap.add_argument('--extension',type=Path,default=EXT,help='Root of the exploratory extension containing legacy_reference and results')
    ap.add_argument('--out',type=Path,default=OUT)
    ap.add_argument('--private-out',type=Path,default=PRIVATE)
    args=ap.parse_args()
    SOURCE=args.source.resolve();EXT=args.extension.resolve();OUT=args.out.resolve();PRIVATE=args.private_out.resolve()
    run()
