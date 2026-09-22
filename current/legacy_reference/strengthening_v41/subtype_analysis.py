#!/usr/bin/env python3
"""Reaggregation and annotation sensitivity, using downloaded original files."""
from pathlib import Path
import argparse, json, re
import h5py, numpy as np, pandas as pd
from scipy import sparse,stats
from strengthen_fc import bh,compare,fit

def run(root, source=None, reuse_checkpoint=False):
    out=root/'results'; private=out/'private'
    allmeta=pd.read_csv(root/'metadata_all.csv.gz',usecols=['barcode','sample','cluster','level2'])
    ic=pd.read_csv(root/'metadata_ic.csv.gz',usecols=['barcode','ic_cluster'])
    merged=allmeta.merge(ic,on='barcode',how='left',validate='one_to_one')
    broad=merged[merged.cluster.isin(['Macro1','Macro2'])].copy()
    broad['ic_cluster']=broad.ic_cluster.fillna('Unassigned')
    counts=pd.crosstab(broad['sample'],broad.ic_cluster)
    d=pd.read_csv(private/'fc_donor_analysis.csv').set_index('sample');counts=counts.loc[d.index]
    fractions=counts.div(counts.sum(axis=1),axis=0)
    rows=[]
    for c in fractions:
        r=compare(fractions[c],d)
        rows.append(dict(refined_annotation=c,**r,total_nuclei=int(counts[c].sum())))
    tab=pd.DataFrame(rows);tab['q_all_refined_categories']=bh(tab.p_mannwhitney)
    tab.to_csv(out/'macrophage_annotation_composition.csv',index=False)
    np.log(counts+.5).pipe(lambda x:x.sub(x.mean(axis=1),axis=0)).to_csv(private/'composition_clr.csv')
    clr=np.log(counts+.5).to_numpy();clr-=clr.mean(axis=1,keepdims=True);clr-=clr.mean(axis=0,keepdims=True)
    u,s,v=np.linalg.svd(clr,full_matrices=False);pc=u[:,:2]*s[:2];d['composition_pc1']=pc[:,0];d['composition_pc2']=pc[:,1]
    r=fit(d.fc_score,d,('age','center','composition_pc1','composition_pc2'))
    comp=dict(model='Age, center and two composition principal components',**r,pc1_variance_fraction=float(s[0]**2/(s**2).sum()),pc2_variance_fraction=float(s[1]**2/(s**2).sum()),n_composition_categories=counts.shape[1])
    pd.DataFrame([comp]).to_csv(out/'fc_composition_adjustment.csv',index=False)
    d.join(fractions.add_prefix('fraction_')).to_csv(private/'fc_composition_donor.csv')
    fc=['FCGR1A','FCGR2A','FCGR2B','FCGR3A','FCGR3B','FCGRT','FCER1G','TYROBP','SYK','LYN','HCK']
    if not reuse_checkpoint:
        allrows=[]; coverage=[]
        for p in sorted((root/'nerve_h5').glob('*.h5')):
            sample=p.name.split('_')[1]
            if sample not in d.index: continue
            ann=merged[merged['sample'].eq(sample)].copy();ann['raw_barcode']=ann.barcode.str.slice(len(sample)+1);ann=ann.set_index('raw_barcode')
            with h5py.File(p) as f:
                m=f['matrix'];genes=[g.decode() for g in m['features/name'][:]];bcs=[x.decode() for x in m['barcodes'][:]]
                X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
                a=ann.reindex(bcs);gt={g:i for i,g in enumerate(genes)}
                Y=X[[gt[g] for g in fc],:]; libpercell=np.asarray(X.sum(axis=0)).ravel()
                masks={'Original macrophages':a.cluster.isin(['Macro1','Macro2']),'Original Macro1':a.cluster.eq('Macro1'),'Original Macro2':a.cluster.eq('Macro2'),'Refined macrophages':a.ic_cluster.str.fullmatch(r'Macro\d+',na=False)}
                masks.update({f'Refined Macro{i}':a.ic_cluster.eq(f'Macro{i}') for i in range(1,19)})
                for label,mask in masks.items():
                    ix=np.flatnonzero(mask.to_numpy());n=len(ix);lib=float(libpercell[ix].sum())
                    coverage.append(dict(sample=sample,disease=d.loc[sample,'disease'],annotation=label,n_nuclei=n,library_size=lib))
                    if n<20 or lib<=0: continue
                    total=np.asarray(Y[:,ix].sum(axis=1)).ravel()
                    for g,c in zip(fc,total): allrows.append(dict(sample=sample,disease=d.loc[sample,'disease'],annotation=label,gene=g,counts=float(c),library_size=lib,n_nuclei=n,log2_cpm_plus_0_5=float(np.log2(c/lib*1e6+.5))))
            print('Aggregated',sample,flush=True)
        ex=pd.DataFrame(allrows);ex.to_csv(private/'fc_refined_pseudobulk.csv',index=False)
        pd.DataFrame(coverage).to_csv(private/'fc_annotation_coverage.csv',index=False)
    else:
        ex=pd.read_csv(private/'fc_refined_pseudobulk.csv')
        coverage=pd.read_csv(private/'fc_annotation_coverage.csv').to_dict('records')
    old=pd.read_csv((source or root.parent/'release/gbs-cidp-blood-nerve-barrier/results')/'tables/gse285983_sample_targeted_expression.csv')
    old=old[old.cell_group.eq('Macrophage') & old.gene.isin(fc) & old['sample'].isin(d.index)]
    check=ex[ex.annotation.eq('Original macrophages')].merge(old,on=['sample','gene'],validate='one_to_one',suffixes=('_new','_old'))
    err=float(np.max(np.abs(check.log2_cpm_plus_0_5_new-check.log2_cpm_plus_0_5_old)))
    assert len(check)==20*11 and err<1e-10,(len(check),err)
    tests=[];eligible=[];score_rows=[]
    for label,part in ex.groupby('annotation',sort=False):
        mat=part.pivot(index='sample',columns='gene',values='log2_cpm_plus_0_5');dd=d.loc[mat.index]
        nc=int(dd.disease.eq('CIDP').sum());nr=int(dd.disease.eq('CIAP').sum())
        eligible.append(dict(annotation=label,n_case=nc,n_reference=nr,eligible=nc>=3 and nr>=3))
        if nc<3 or nr<3:continue
        z=(mat-mat.mean())/mat.std(ddof=1).replace(0,np.nan);scores=z.mean(axis=1)
        tests.append(dict(annotation=label,**compare(scores,dd)))
        score_rows.extend(dict(sample=i,annotation=label,score=float(scores.loc[i])) for i in scores.index)
    tab2=pd.DataFrame(tests); subtype=tab2.annotation.str.fullmatch(r'Refined Macro\d+')
    tab2['q_eligible_refined_subtypes']=np.nan;tab2.loc[subtype,'q_eligible_refined_subtypes']=bh(tab2.loc[subtype,'p_mannwhitney'])
    gate=~subtype;tab2['q_eligible_annotation_definitions']=np.nan;tab2.loc[gate,'q_eligible_annotation_definitions']=bh(tab2.loc[gate,'p_mannwhitney'])
    tab2.to_csv(out/'fc_annotation_sensitivity.csv',index=False)
    seen={x['annotation'] for x in eligible}
    for label in pd.DataFrame(coverage).annotation.unique():
        if label not in seen: eligible.append(dict(annotation=label,n_case=0,n_reference=0,eligible=False))
    pd.DataFrame(eligible).to_csv(out/'fc_subtype_eligibility.csv',index=False)
    pd.DataFrame(score_rows).to_csv(private/'fc_subtype_scores.csv',index=False)
    # Spatial gene panel coverage is a measurement audit, not a disease test.
    spatial=[];gene_sets=[]
    for p in sorted((root/'spatial').glob('*.h5')):
        with h5py.File(p) as f:
            m=f['matrix'];gn=[g.decode() for g in m['features/name'][:]]
            bio=[g for g in gn if not g.startswith(('NegControl','UnassignedCodeword'))];gene_sets.append(set(bio))
            sample=p.name.split('_')[1]; dis=allmeta.loc[allmeta['sample'].eq(sample),'level2'].iloc[0]
            spatial.append(dict(sample=sample,disease=dis,n_segmented_cells=int(m['shape'][1]),n_biological_genes=len(bio),n_fc_panel_genes=len(set(bio)&set(fc+['FCGR2C']))))
    assert all(g==gene_sets[0] for g in gene_sets)
    genes_report=pd.DataFrame([dict(gene=g,measured_in_spatial=g in gene_sets[0]) for g in fc+['FCGR2C']]);genes_report.to_csv(out/'spatial_fc_coverage.csv',index=False)
    pd.DataFrame(spatial).to_csv(private/'spatial_sample_coverage.csv',index=False)
    pd.DataFrame(spatial).groupby('disease').agg(n_donors=('sample','size'),n_segmented_cells=('n_segmented_cells','sum'),n_biological_genes=('n_biological_genes','min'),n_fc_genes=('n_fc_panel_genes','min')).reset_index().to_csv(out/'spatial_design_summary.csv',index=False)
    summary={'raw_reproduction_max_error':err,'raw_reproduced_values':len(check),'refined_category_tests':len(tab),'eligible_refined_subtype_tests':int(subtype.sum()),'composition_model':comp,'spatial_fc_coverage':0,'spatial_biological_panel':len(gene_sets[0]),'annotation_tests':tab2.replace({np.nan:None}).to_dict('records')}
    (out/'subtype_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2));print(tab.nsmallest(5,'p_mannwhitney')[['refined_annotation','delta','p_mannwhitney','q_all_refined_categories']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--original-results',type=Path);p.add_argument('--reuse-checkpoint',action='store_true');a=p.parse_args();run(a.root,a.original_results,a.reuse_checkpoint)
