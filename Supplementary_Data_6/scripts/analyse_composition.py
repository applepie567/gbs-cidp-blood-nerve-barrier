"""Donor- and centre-held-out evaluation, followed by a gated external analysis."""
from pathlib import Path
import sys, json, re
import numpy as np
import pandas as pd
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'input/prior/External_nerve_20260926'
sys.path.insert(0,str(PRIOR))
from analyse_external_nerve import GENES, hc3
OUT=ROOT/'results'
SEX={'XIST','RPS4Y1','RPS4Y2','DDX3Y','EIF1AY','KDM5D','UTY','ZFY','USP9Y','TMSB4Y','NLGN4Y','TXLNGY'}

def corr(x,y):
    if len(x)<3 or np.std(x)==0 or np.std(y)==0: return float('nan')
    return float(stats.spearmanr(x,y).statistic)

def run():
    z=np.load(OUT/'atlas_donor_lineage_counts.npz')
    counts=z['counts']; genes=z['genes'].astype(str)
    rm=pd.read_csv(OUT/'atlas_donor_lineage_metadata.csv')
    assert counts.shape[0]==len(rm)
    libs=rm.library_counts.to_numpy()
    cpm=np.divide(counts,libs[:,None],out=np.zeros_like(counts),where=libs[:,None]>0)*1e6
    probes=pd.read_csv(PRIOR/'inputs/GSE213455_series_matrix.txt.gz',sep='\t',comment='!',index_col=0)
    probes.index=probes.index.astype(str)
    annotation_path=PRIOR/'inputs/GPL13369_annotation.tsv'
    if not annotation_path.exists(): annotation_path=annotation_path.with_suffix('.tsv.gz')
    ann=pd.read_csv(annotation_path,sep='\t',usecols=['ID','SYMBOL','CHROMOSOME'],dtype=str).fillna('')
    ann=ann[ann.ID.isin(probes.index)&ann.SYMBOL.str.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*')].copy()
    ambiguous=ann.groupby('ID').SYMBOL.nunique()
    ann=ann[ann.ID.isin(ambiguous[ambiguous==1].index)].drop_duplicates('ID')
    ann.to_csv(OUT/'external_unambiguous_probe_mapping.csv',index=False)
    common=sorted(set(genes)&set(ann.SYMBOL)-set(GENES))
    ix=np.array([np.flatnonzero(genes==g)[0] for g in common])
    pd.DataFrame({'gene':common}).to_csv(OUT/'rank_universe.csv',index=False)
    candidates=np.array([g not in SEX and not re.match(r'^(MT-|RPL|RPS)',g) for g in common])
    lineages=rm.lineage.unique().tolist(); donors=sorted(rm.sample_id.unique())
    dm=rm[['sample_id','diagnosis','age','sex','centre']].drop_duplicates().set_index('sample_id').loc[donors]
    whole=[]; targets=[]
    for sample in donors:
        take=rm.sample_id==sample
        co=counts[take].sum(axis=0)
        whole.append(co[ix])
        ma=rm[take&(rm.lineage=='Macrophages')].iloc[0]
        targets.append(dict(sample_id=sample,nucleus_fraction=ma.n_nuclei/rm.loc[take,'n_nuclei'].sum(),
                            captured_count_fraction=ma.library_counts/rm.loc[take,'library_counts'].sum()))
    whole=np.array(whole)
    ranks=stats.rankdata(whole,axis=1,method='average')/len(common)
    tg=pd.DataFrame(targets).set_index('sample_id').join(dm)

    def select(excluded):
        means=[]; eligible=rm.reference_eligible.astype(bool)&~rm.sample_id.isin(excluded)
        for lin in lineages:
            e=eligible&(rm.lineage==lin)
            # Each observed lineage is retained; no missing profile is imputed.
            means.append(cpm[e][:,ix].mean(axis=0) if e.sum() else np.zeros(len(ix)))
        means=np.array(means); mi=lineages.index('Macrophages')
        em=eligible&(rm.lineage=='Macrophages')
        detected=(cpm[em][:,ix]>1).mean(axis=0)
        enrichment=(means[mi]+1)/(np.delete(means,mi,axis=0).max(axis=0)+1)
        ok=candidates&(means[mi]>=10)&(detected>=.5)&(enrichment>=4)
        k=np.flatnonzero(ok)
        k=k[np.lexsort((np.array(common)[k],-enrichment[k]))][:20]
        table=pd.DataFrame(dict(gene=np.array(common)[k],macrophage_mean_CPM=means[mi,k],
                                enrichment=enrichment[k],donor_detection_fraction=detected[k]))
        return k,table

    selections=[];pred=[]
    for i,sample in enumerate(donors):
        k,table=select([sample]);table['held_out']=sample;table['scheme']='donor'
        selections.append(table)
        pred.append(dict(sample_id=sample,index=float(ranks[i,k].mean()) if len(k)>=10 else np.nan,n_markers=len(k)))
    loo=pd.DataFrame(pred).set_index('sample_id').join(tg)
    loo.to_csv(OUT/'donor_holdout_predictions.csv')
    centre_rows=[]
    for centre in sorted(dm.centre.unique()):
        excluded=dm.index[dm.centre==centre].tolist()
        k,table=select(excluded);table['held_out']=centre;table['scheme']='centre';selections.append(table)
        for sample in excluded:
            i=donors.index(sample)
            centre_rows.append(dict(sample_id=sample,index=float(ranks[i,k].mean()) if len(k)>=10 else np.nan,n_markers=len(k)))
    ch=pd.DataFrame(centre_rows).set_index('sample_id').join(tg)
    ch.to_csv(OUT/'centre_holdout_predictions.csv')
    k,table=select([]);table['held_out']='none';table['scheme']='full_reference';selections.append(table)
    pd.concat(selections,ignore_index=True).to_csv(OUT/'all_marker_selections.csv',index=False)
    table.to_csv(OUT/'frozen_macrophage_markers.csv',index=False)
    metrics=[]
    for scheme,frame in [('donor',loo),('centre',ch)]:
        for centre in ['All']+sorted(dm.centre.unique()):
            s=frame if centre=='All' else frame[frame.centre==centre]
            metrics.append(dict(scheme=scheme,centre=centre,n_donors=len(s),minimum_markers=int(s.n_markers.min()),
                                rho_nucleus_fraction=corr(s['index'],s.nucleus_fraction),
                                rho_captured_count_fraction=corr(s['index'],s.captured_count_fraction)))
    mt=pd.DataFrame(metrics);mt.to_csv(OUT/'feasibility_metrics.csv',index=False)
    primary=mt[(mt.scheme=='donor')&(mt.centre=='All')].iloc[0]
    centre_valid=mt[mt.centre!='All']
    gate=bool(loo['index'].notna().all() and ch['index'].notna().all() and len(k)>=10 and
              primary.rho_captured_count_fraction>=.70 and primary.rho_nucleus_fraction>=.60 and
              (centre_valid[['rho_nucleus_fraction','rho_captured_count_fraction']]>0).all().all())
    summary=dict(n_donors=len(donors),n_nuclei=int(rm.n_nuclei.sum()),n_genes_in_rank_universe=len(common),
                 n_frozen_markers=len(k),excluded_outcome_genes=GENES,gate_passed=gate,
                 primary_rho_nuclei=float(primary.rho_nucleus_fraction),primary_rho_counts=float(primary.rho_captured_count_fraction),
                 interpretation='Marker expression index; not an absolute cell fraction. Atlas benchmarking does not validate cross-platform accuracy.')
    if gate:
        # Only now inspect the existing external Fc outcome against the frozen index.
        em=probes.loc[ann.ID].copy();em['gene']=ann.set_index('ID').loc[em.index,'SYMBOL'].to_numpy()
        external=em.groupby('gene').median().loc[common].T
        external_mean=em.groupby('gene').mean().loc[common].T
        fixed=pd.read_csv(PRIOR/'results/Fc_donor_scores.csv',index_col=0)
        external=external.loc[fixed.index];external_mean=external_mean.loc[fixed.index]
        er=stats.rankdata(external.to_numpy(),axis=1,method='average')/len(common)
        ermean=stats.rankdata(external_mean.to_numpy(),axis=1,method='average')/len(common)
        fixed['macrophage_index']=er[:,k].mean(axis=1)
        fixed['macrophage_index_mean_probe']=ermean[:,k].mean(axis=1)
        fixed.to_csv(OUT/'external_donor_scores.csv')
        external.to_csv(OUT/'external_gene_expression_median.csv.gz')
        external_mean.to_csv(OUT/'external_gene_expression_mean.csv.gz')
        d=(fixed.diagnosis=='CIDP').to_numpy().astype(int);y=fixed.Fc_score.to_numpy()
        cov=stats.zscore(fixed.macrophage_index.to_numpy(),ddof=1)
        models=[]
        for name,yy,cc in [('Unadjusted',y,[]),('Macrophage index adjusted',y,[cov]),
                           ('Mean probe sensitivity',fixed.Fc_score_mean_probe.to_numpy(),[stats.zscore(fixed.macrophage_index_mean_probe,ddof=1)])]:
            rr=hc3(yy,d,cc);rr['model']=name;models.append(rr)
        pd.DataFrame(models).to_csv(OUT/'external_conditional_models.csv',index=False)
        omissions=[]
        for sample in fixed.index:
            sub=fixed.drop(sample);dd=(sub.diagnosis=='CIDP').to_numpy().astype(int)
            rr=hc3(sub.Fc_score.to_numpy(),dd,[stats.zscore(sub.macrophage_index,ddof=1)])
            rr['omitted_sample']=sample;omissions.append(rr)
        pd.DataFrame(omissions).to_csv(OUT/'external_leave_one_out.csv',index=False)
        X=np.column_stack([np.ones(len(d)),d,cov]);hat=np.diag(X@np.linalg.inv(X.T@X)@X.T)
        fixed['conditional_leverage']=hat
        fixed.to_csv(OUT/'external_donor_scores.csv')
        summary['external']=dict(n_CIDP=4,n_VN=9,rho_index_Fc=corr(fixed.macrophage_index,fixed.Fc_score),
             max_leverage=float(hat.max()),models=models,
             omitted_coefficient_min=float(min(x['difference'] for x in omissions)),omitted_coefficient_max=float(max(x['difference'] for x in omissions)),
             marker_index_CIDP_range=[float(fixed.loc[fixed.diagnosis=='CIDP','macrophage_index'].min()),float(fixed.loc[fixed.diagnosis=='CIDP','macrophage_index'].max())],
             marker_index_VN_range=[float(fixed.loc[fixed.diagnosis=='VN','macrophage_index'].min()),float(fixed.loc[fixed.diagnosis=='VN','macrophage_index'].max())])
    else:
        summary['external']='Not fitted because the atlas feasibility gate failed.'
    (OUT/'composition_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False))
    print(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False),flush=True)
    print(mt.to_string(index=False),flush=True)

if __name__=='__main__':run()
