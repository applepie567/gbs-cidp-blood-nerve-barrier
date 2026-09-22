"""Exploratory partition of nerve programs, with donor-level inference.

Run from any directory with explicit --raw, --legacy and --private paths.
The input manifest and analysis_plan.md document provenance and selection.
"""
import argparse
import gc
import hashlib
import json
import platform
import re
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats
from read_h5 import H5Reader
from stats_utils import bh, hc3

PANELS = {
    'Endothelial junctions': ['CLDN5','OCLN','TJP1','CDH5'],
    'Endothelial transport': ['ABCB1','SLC1A1','MFSD2A'],
    'Inflammatory adhesion': ['ICAM1','VCAM1','SELE','SELP'],
    'Myelin maintenance': ['MPZ','MBP','PRX','EGR2','PMP22','MAG'],
    'Injury response': ['NGFR','ATF3','JUN','FOS','GDNF','RUNX2'],
}
FC = ['FCGR1A','FCGR2A','FCGR2B','FCGR2C','FCGR3A','FCGR3B','FCGRT','FCER1G','TYROBP','SYK','LYN','HCK']
GROUPS = {'BNB_EC':['ven_capEC2'], 'Myelinating_SC':['mySC'],
          'Nonmyelinating_SC':['nmSC'], 'Repair_damage_SC':['repairSC','damageSC'],
          'Macrophage':['Macro1','Macro2']}
TARGETS = {'BNB_EC': list(PANELS)[:3], 'Myelinating_SC': list(PANELS)[3:],
           'Nonmyelinating_SC': list(PANELS)[3:], 'Repair_damage_SC': list(PANELS)[3:]}


def run(raw, legacy, out, private):
    out.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    cols = ['sample','barcode','cluster','level2','sex','age','center','incat']
    md = pd.read_csv(raw/'metadata_all.csv.gz', usecols=cols, dtype=str)
    md = md[md.level2.isin(['CIDP','CIAP'])].copy()
    donor = md.groupby('sample')[['level2','sex','age','center','incat']].nunique(dropna=False)
    assert donor.to_numpy().max() == 1, 'Donor metadata not constant'
    donors = md.drop_duplicates('sample').set_index('sample')
    assert donors.level2.value_counts().to_dict() == {'CIAP':11, 'CIDP':9}
    genes = sorted(set(sum(PANELS.values(), []) + FC + ['SOX10','S100B']))
    long, elig, reconstruction = [], [], []
    backends = set()
    for path in sorted((raw/'nerve_h5').glob('*.h5')):
        sample = re.search(r'_(S\d+)_', path.name).group(1)
        if sample not in donors.index:
            continue
        dm = md[md['sample'].eq(sample)].copy()
        assert dm.barcode.str.startswith(sample+'_').all()
        dm['raw_barcode'] = dm.barcode.str.slice(len(sample)+1)
        assert not dm.raw_barcode.duplicated().any()
        dm = dm.set_index('raw_barcode')
        meta = donors.loc[sample]
        with H5Reader(path) as f:
            backends.add(f.backend)
            names = f.array('matrix/features/name')
            barcodes = f.array('matrix/barcodes')
            shape = tuple(int(x) for x in f.array('matrix/shape'))
            indptr = f.array('matrix/indptr')
            indices = f.array('matrix/indices')
            data = f.array('matrix/data')
        assert len(names)==shape[0] and len(barcodes)==shape[1]
        assert len(set(barcodes))==len(barcodes)
        assert len(data)==len(indices)==int(indptr[-1])
        assert (data>=0).all() and np.isfinite(data).all()
        assert indices.max()<shape[0] and np.all(np.diff(indptr)>=0)
        assert set(dm.index).issubset(set(barcodes)), 'Unmatched annotated nuclei'
        mat = sparse.csc_matrix((data,indices,indptr), shape=shape, copy=False)
        lookup = {g: np.flatnonzero(names==g) for g in genes}
        present = [g for g in genes if len(lookup[g])]
        assert all(len(lookup[g])==1 for g in present), 'Duplicate target gene'
        rows = np.array([lookup[g][0] for g in present])
        targeted = mat[rows,:]
        annotation = dm.reindex(barcodes)['cluster']
        for group, labels in GROUPS.items():
            take = np.flatnonzero(annotation.isin(labels).to_numpy())
            n = len(take)
            fraction = float(annotation.iloc[take].eq('repairSC').mean()) if group=='Repair_damage_SC' and n else np.nan
            elig.append(dict(sample=sample,disease=meta.level2,cell_group=group,n_nuclei=n,
                             eligible=n>=20,repair_fraction=fraction))
            if n<20:
                continue
            lib = float(mat[:,take].sum())
            total = np.asarray(targeted[:,take].sum(axis=1)).ravel()
            assert lib>0 and np.sum(total)<=lib
            for g,count in zip(present,total):
                long.append(dict(sample=sample,disease=meta.level2,sex=meta.sex,age=float(meta.age),
                                 center=str(meta.center),incat=float(meta.incat) if pd.notna(meta.incat) else np.nan,
                                 cell_group=group,n_nuclei=n,library_size=lib,gene=g,count=float(count),
                                 log2_cpm_plus_0_5=float(np.log2(count/lib*1e6+.5)),repair_fraction=fraction))
        reconstruction.append(dict(sample=sample,annotated_nuclei=len(dm),matrix_columns=shape[1],
                                   matrix_genes=shape[0],target_genes_present=len(present),nnz=len(data)))
        print('Aggregated',sample,meta.level2,'annotated nuclei',len(dm),flush=True)
        del mat,targeted,data,indices,indptr
        gc.collect()
    e = pd.DataFrame(long)
    e.to_csv(private/'nerve_targeted_expression.csv',index=False)
    pd.DataFrame(elig).to_csv(private/'nerve_donor_eligibility.csv',index=False)
    pd.DataFrame(reconstruction).to_csv(private/'nerve_reconstruction.csv',index=False)
    return compute(e, pd.DataFrame(elig), backends, legacy, out, private)


def compute(e, elig, backends, legacy, out, private):
    coverage, score_rows, results, adjusted, composition, gene_summaries = [],[],[],[],[],[]
    for group, programs in TARGETS.items():
        x = e[e.cell_group.eq(group)]
        d = x.drop_duplicates('sample').set_index('sample').sort_index()
        m = x.pivot(index='sample',columns='gene',values='log2_cpm_plus_0_5').loc[d.index]
        for prog in programs:
            usable=[]
            for g in PANELS[prog]:
                exists=g in m.columns
                var=float(m[g].std(ddof=1)) if exists else np.nan
                use=exists and np.isfinite(var) and var>0
                coverage.append(dict(cell_group=group,program=prog,gene=g,present=exists,
                                     sd_across_donors=var,in_score=use,
                                     n_donors_detected=int((m[g]>-1+1e-10).sum()) if exists else 0,
                                     n_eligible_donors=len(d)))
                if use: usable.append(g)
                if exists:
                    for disease in ['CIDP','CIAP']:
                        vals=m.loc[d.index[d.disease.eq(disease)],g]
                        gene_summaries.append(dict(cell_group=group,program=prog,gene=g,disease=disease,
                            n_donors=len(vals),mean_log2_cpm=float(vals.mean()),sd_log2_cpm=float(vals.std(ddof=1)),
                            detected_donors=int((vals>-1+1e-10).sum())))
            assert len(usable)>=2, (group,prog,usable)
            z=(m[usable]-m[usable].mean())/m[usable].std(ddof=1)
            y=z.mean(axis=1)
            for sample,value in y.items():
                score_rows.append(dict(sample=sample,disease=d.loc[sample,'disease'],cell_group=group,program=prog,
                                       score=value,incat=d.loc[sample,'incat'],n_nuclei=int(d.loc[sample,'n_nuclei'])))
            case=d.disease.eq('CIDP')
            r=dict(cell_group=group,program=prog,genes_planned=len(PANELS[prog]),genes_used=len(usable),
                   n_case=int(case.sum()),n_reference=int((~case).sum()),mean_case=float(y[case].mean()),
                   mean_reference=float(y[~case].mean()),p_mannwhitney=float(stats.mannwhitneyu(y[case],y[~case],alternative='two-sided').pvalue),
                   **hc3(y,d))
            results.append(r)
            adjusted.append(dict(cell_group=group,program=prog,model='Age and centre',**hc3(y,d,['age','center'])))
            if group=='Repair_damage_SC':
                composition.append(dict(cell_group=group,program=prog,model='Age centre and repairSC fraction',
                                        **hc3(y,d,['age','center','repair_fraction'])))
    res=pd.DataFrame(results)
    assert len(res)==9
    res['q_bh_9']=bh(res.p_mannwhitney)
    old=pd.read_csv(legacy/'strengthening_v41/results/module_global_multiplicity.csv')
    assert len(old)==80
    family=bh(np.r_[old.p_value.to_numpy(),res.p_mannwhitney.to_numpy()])
    old['q_bh_expanded_89']=family[:80]
    res['q_bh_89']=family[80:]
    adj=pd.DataFrame(adjusted)
    adj['q_bh_9_hc3']=bh(adj.p_hc3)
    comp=pd.DataFrame(composition)
    comp['q_bh_2_hc3']=bh(comp.p_hc3)
    scores=pd.DataFrame(score_rows)
    scores.to_csv(private/'nerve_program_scores.csv',index=False)
    corr=[]
    for prog in list(PANELS)[3:]:
        d=scores.query('cell_group=="Repair_damage_SC" and disease=="CIDP"').query('program==@prog').dropna(subset=['incat','score'])
        rho,p=stats.spearmanr(d.score,d.incat)
        corr.append(dict(cell_group='Repair_damage_SC',program=prog,n_donors=len(d),rho=rho,p_value=p))
    cor=pd.DataFrame(corr)
    cor['q_bh_2']=bh(cor.p_value)
    oldc=pd.read_csv(legacy/'analysis_update/source_data/CIDP_clinical.csv')
    assert len(oldc)==15
    qc=bh(np.r_[oldc.p_value.to_numpy(),cor.p_value.to_numpy()])
    oldc['q_bh_expanded_17']=qc[:15]
    cor['q_bh_17']=qc[15:]
    cor=cor.rename(columns={'p_value':'p_value_asymptotic','q_bh_2':'q_asymptotic_BH_2','q_bh_17':'q_asymptotic_BH_17'})
    # Processing check against an archived raw-matrix analysis with the same inputs.
    fm=e.query('cell_group=="Macrophage"').pivot(index='sample',columns='gene',values='log2_cpm_plus_0_5')
    fd=e.query('cell_group=="Macrophage"').drop_duplicates('sample').set_index('sample').loc[fm.index]
    measured=[g for g in FC if g in fm.columns]
    fs=((fm[measured]-fm[measured].mean())/fm[measured].std(ddof=1).replace(0,np.nan)).mean(axis=1)
    fdelta=float(fs[fd.disease.eq('CIDP')].mean()-fs[fd.disease.eq('CIAP')].mean())
    fp=float(stats.mannwhitneyu(fs[fd.disease.eq('CIDP')],fs[fd.disease.eq('CIAP')],alternative='two-sided').pvalue)
    prior=pd.read_csv(legacy/'strengthening_v41/results/fc_sensitivity.csv')
    p20=prior[prior.analysis.eq('Standardisation in CIDP and CIAP only')].iloc[0]
    np.testing.assert_allclose([fdelta,fp],[p20.delta,p20.p_primary],atol=1e-10)
    pg=pd.read_csv(legacy/'strengthening_v41/results/fc_component_genes.csv').set_index('gene')
    discrepancies=[]
    for g in measured:
        delta=float(fm.loc[fd.disease.eq('CIDP'),g].mean()-fm.loc[fd.disease.eq('CIAP'),g].mean())
        discrepancies.append(abs(delta-pg.loc[g,'delta']))
    assert max(discrepancies)<1e-10
    for name,df in [('nerve_program_contrasts',res),('nerve_program_adjusted',adj),
                    ('schwann_composition_sensitivity',comp),('nerve_gene_coverage',pd.DataFrame(coverage)),
                    ('nerve_gene_summary',pd.DataFrame(gene_summaries)),('nerve_clinical_correlations',cor),
                    ('original_nerve_multiplicity_89',old),('original_clinical_asymptotic_17',oldc)]:
        df.to_csv(out/(name+'.csv'),index=False)
    el=pd.DataFrame(elig)
    el.groupby(['disease','cell_group']).agg(donors=('sample','nunique'),eligible_donors=('eligible','sum'),
        nuclei_min=('n_nuclei','min'),nuclei_max=('n_nuclei','max')).reset_index().to_csv(out/'nerve_eligibility_summary.csv',index=False)
    summary=dict(exploratory=True,normalisation='log2(CPM + 0.5), gene z scores across eligible CIDP and CIAP donors',
                 case_donors=9,comparator_donors=11,new_disease_tests=9,expanded_disease_tests=89,
                 new_clinical_tests=2,expanded_clinical_tests=17,fc_reconstruction_delta=fdelta,
                 fc_reconstruction_p=fp,max_fc_gene_delta_error=max(discrepancies),
                 hdf5_backends=sorted(backends),versions={'python':platform.python_version(),'numpy':np.__version__,
                 'pandas':pd.__version__,'scipy':scipy.__version__})
    (out/'nerve_summary.json').write_text(json.dumps(summary,indent=2))
    print(res[['cell_group','program','delta','ci_low','ci_high','p_mannwhitney','q_bh_9','q_bh_89']].to_string(index=False))
    print(cor.to_string(index=False))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['raw','legacy','out','private']:
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--reuse-reconstruction', action='store_true')
    a=p.parse_args()
    if a.reuse_reconstruction:
        with H5Reader(next((a.raw/'nerve_h5').glob('*.h5'))) as f:
            backends={f.backend}
        compute(pd.read_csv(a.private/'nerve_targeted_expression.csv'),
                pd.read_csv(a.private/'nerve_donor_eligibility.csv'),backends,a.legacy,a.out,a.private)
    else:
        run(a.raw,a.legacy,a.out,a.private)
