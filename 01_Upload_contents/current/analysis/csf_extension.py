"""Reproducible synthesis of printed, adjusted CSF coefficients.

Source: Kmezic et al. Front Immunol. 2023;14:1241199.
https://doi.org/10.3389/fimmu.2023.1241199
All inputs below are transcribed from the source's printed Tables 3, 4B, 5, 6.
The individual-level NPX values are not available here and are not reconstructed.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from stats_utils import bh, meta_two

URL='https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1241199/full'
# protein, discovery beta, SE, P, FDR, replication beta, SE, P, FDR
T4B=[
 ['SELE',1.79,.377,.0004,.0086,2.33,.766,.0188,.0446],
 ['IL2RA',.741,.234,.0075,.0466,2.47,.313,9.94e-5,.0019],
 ['CCL3',.836,.270,.0085,.0466,.684,.534,.241,.294],
 ['CD4',.359,.131,.0169,.0743,None,None,None,None],
 ['MMP7',1.2,.470,.0245,.0786,None,None,None,None],
 ['CR2',.692,.273,.025,.0786,1.92,.624,.0178,.0446],
 ['CD1C',.561,.252,.0446,.123,1.04,.299,.0102,.0388],
 ['THBD',.656,.317,.059,.144,1.92,.477,.0051,.0238],
 ['NRP1',.0965,.0485,.0683,.15,.288,.119,.0459,.0793],
 ['CD38',.635,.352,.0943,.173,1.18,.447,.0389,.0738],
 ['IL3RA',.183,.0982,.0867,.173,None,None,None,None],
 ['IL6',.37,.514,.485,.547,1.89,.427,.0031,.0194],
 ['ITGAM',None,None,None,None,1.15,.351,.0136,.0429],
 ['CD5',-.0882,.126,.497,.547,1.51,.514,.0218,.046],
]
# Every CSF IL8 contrast printed in Table 3, fully adjusted model.
T3=[
 ['GBS versus HC','Discovery',1.66,.415,.0015,.0163,6,20],
 ['GBS versus HC','Replication',2.62,.462,.0008,.0071,6,6],
 ['GBS versus CIDP','Replication',1.71,.419,.0046,.0291,6,6],
 ['GBS versus NIP','Replication',2.96,.837,.0167,.0444,6,4],
 ['GBS versus ALS','Discovery',1.61,.342,.0022,.0159,6,6],
 ['CIDP versus HC','Replication',1.13,.181,.0004,.0081,6,6],
 ['CIDP versus NIP','Replication',2.08,.930,.0759,.406,6,4],
]
# protein, fully adjusted beta, SE, source P, source FDR
T5=[
 ['ITGAM',1.19,.26,.0026,.0291], ['IL2RA',1.58,.427,.0077,.0364],
 ['IL6',1.33,.315,.0039,.0291], ['NRP1',.387,.114,.0116,.0441],
 ['SELE',1.91,.785,.0455,.108], ['THBD',1.09,.481,.0577,.11],
 ['CD40',.519,.183,.0254,.0804], ['CR2',1.56,.669,.0527,.11],
 ['CD28',.557,.26,.0693,.12],
]
T6=[
 ['MMP3',.57,.238,.0477,.181], ['THBD',1.04,.397,.0341,.181],
 ['IL2RA',1.16,.458,.0387,.181], ['CR2',.893,.493,.113,.315],
 ['CD1C',.765,.303,.0396,.181],
]


def run(out, source, legacy):
    out.mkdir(parents=True,exist_ok=True)
    source.mkdir(parents=True,exist_ok=True)
    rows=[]
    for protein,*values in T4B:
        for k,cohort in enumerate(['Discovery','Replication']):
            beta,se,p,q=values[k*4:k*4+4]
            rows.append(dict(protein=protein,contrast='GBS versus HC',cohort=cohort,
                beta=beta,se=se,source_p=p,source_q=q,n_case=6,n_reference=20 if k==0 else 6,
                source_table='4B',availability='reported' if beta is not None else 'NA in source due to call rate below 60%'))
    for contrast,cohort,beta,se,p,q,nc,nr in T3:
        rows.append(dict(protein='IL8',contrast=contrast,cohort=cohort,beta=beta,se=se,source_p=p,
                         source_q=q,n_case=nc,n_reference=nr,source_table='3',availability='reported'))
    for table,contrast,content in [('5','GBS versus CIDP',T5),('6','CIDP versus HC',T6)]:
        for protein,beta,se,p,q in content:
            rows.append(dict(protein=protein,contrast=contrast,cohort='Replication',beta=beta,se=se,
                         source_p=p,source_q=q,n_case=6,n_reference=6,source_table=table,availability='reported'))
    src=pd.DataFrame(rows)
    src['specimen']='CSF'
    src['scale']='log2 NPX'
    src['adjustment']='Age sex and CCL19 sample handling marker'
    src['source_url']=URL
    src['verification_date']='2026-09-21'
    src['normal_ci_low']=src.beta-stats.norm.ppf(.975)*src.se
    src['normal_ci_high']=src.beta+stats.norm.ppf(.975)*src.se
    src['interval_method']='Normal approximation from printed SE. Source P and q retained separately.'
    src.to_csv(source/'CSF_all_reported_adjusted_coefficients.csv',index=False)
    paired=[]
    eligibility=[]
    for protein,d in src.query('contrast=="GBS versus HC"').groupby('protein',sort=False):
        complete=d.dropna(subset=['beta','se'])
        eligible=len(complete)==2
        eligibility.append(dict(protein=protein,reported_cohorts=len(complete),included_in_meta=eligible,
            reason='Both independent coefficients available' if eligible else 'One cohort coefficient unavailable'))
        if not eligible:continue
        dd=d.set_index('cohort')
        r=dict(protein=protein,source_table=','.join(sorted(set(d.source_table))),
               discovery_beta=dd.loc['Discovery','beta'],discovery_se=dd.loc['Discovery','se'],
               replication_beta=dd.loc['Replication','beta'],replication_se=dd.loc['Replication','se'],
               discovery_source_p=dd.loc['Discovery','source_p'],discovery_source_q=dd.loc['Discovery','source_q'],
               replication_source_p=dd.loc['Replication','source_p'],replication_source_q=dd.loc['Replication','source_q'],
               **meta_two(d.beta.to_numpy(),d.se.to_numpy()))
        paired.append(r)
    result=pd.DataFrame(paired)
    assert len(result)==11
    result['q_fixed_BH_11']=bh(result.fixed_p)
    result['q_random_BH_11']=bh(result.random_p)
    result['q_mhk_BH_11']=bh(result.mhk_p)
    result['positive_in_both']=(result.discovery_beta>0)&(result.replication_beta>0)
    result['source_q_below_05_in_both']=(result.discovery_source_q<.05)&(result.replication_source_q<.05)
    result.to_csv(out/'CSF_meta_11_proteins.csv',index=False)
    pd.DataFrame(eligibility).to_csv(out/'CSF_extraction_eligibility.csv',index=False)
    compare=src[src.contrast.isin(['GBS versus CIDP','CIDP versus HC'])].copy()
    assert len(compare)==16
    compare['shared_participants']='Overlapping contrasts in the same replication cohort. Not pooled together.'
    compare.to_csv(out/'CSF_disease_comparisons.csv',index=False)
    # The prior ten inputs and point estimates must remain unchanged.
    old=pd.read_csv(legacy/'analysis_update/source_data/CSF_Kmezic_published_coefficients.csv').set_index('protein')
    common=result.set_index('protein').loc[old.index]
    for col in ['discovery_beta','discovery_se','replication_beta','replication_se']:
        np.testing.assert_allclose(common[col],old[col],rtol=0,atol=1e-12)
    summary=dict(exploratory=True,source_tables=['3','4B','5','6'],unique_gbs_hc_proteins=len(eligibility),
        paired_proteins=11,unpaired_proteins=[r['protein'] for r in eligibility if not r['included_in_meta']],
        positive_in_both=int(result.positive_in_both.sum()),fixed_q_below_05=int((result.q_fixed_BH_11<.05).sum()),
        source_q_below_05_in_both=result.loc[result.source_q_below_05_in_both,'protein'].tolist(),
        random_ci_crosses_zero=result.loc[(result.random_low<=0)&(result.random_high>=0),'protein'].tolist(),
        mhk_ci_crosses_zero=result.loc[(result.mhk_low<=0)&(result.mhk_high>=0),'protein'].tolist(),
        source_individual_models_refit=False,il8=result[result.protein.eq('IL8')].iloc[0].to_dict())
    (out/'CSF_summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:bool(x) if isinstance(x,np.bool_) else float(x)))
    print(result[['protein','fixed_beta','fixed_low','fixed_high','q_fixed_BH_11','random_low','random_high','mhk_low','mhk_high']].to_string(index=False))
    print(json.dumps(summary,indent=2,default=lambda x:bool(x) if isinstance(x,np.bool_) else float(x)))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['out','source','legacy']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    run(a.out,a.source,a.legacy)
