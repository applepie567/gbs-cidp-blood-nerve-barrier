"""Offline checks of figure inputs and reported statistical summaries."""
from pathlib import Path
import argparse, hashlib, itertools, json, math, zipfile
import numpy as np,pandas as pd
from scipy import stats
from stats_utils import bh,meta_two

ROOT=Path(__file__).resolve().parents[1]
def rd(path):return pd.read_csv(ROOT/path)
def close(a,b):np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-10,equal_nan=True)
def meta(g,v):
    k=len(g);w=1/v;m=np.average(g,weights=w);q=np.sum(w*(g-m)**2)
    tau=max(0,(q-k+1)/(w.sum()-(w*w).sum()/w.sum()))
    wr=1/(v+tau);m=np.average(g,weights=wr);se=np.sqrt(max(1,np.sum(wr*(g-m)**2)/(k-1))/wr.sum())
    t=stats.t.ppf(.975,k-1)
    return [m,m-t*se,m+t*se,2*stats.t.sf(abs(m/se),k-1)]

def main():
    report={}
    index=rd('figure_source_data/INDEX.csv')
    for r in index.itertuples():assert len(rd('figure_source_data/'+r.csv))==r.rows,r.csv
    expected={f'Figure {i}' for i in range(1,6)}|{f'Supplementary Figure {i}' for i in range(1,8)}
    assert expected.issubset(set(index.figure));report['figure_coverage']=12;report['source_tables']=len(index)
    report['image_hashes_matched']=0
    for r in json.loads((ROOT/'metadata/embedded_figure_manifest.json').read_text()):
        data=(ROOT/'figures'/f"{r['figure']}.png").read_bytes()
        assert hashlib.sha256(data).hexdigest()==r['sha256'],r['figure'];report['image_hashes_matched']+=1
    b=rd('figure_source_data/F2A_Blood_cohort.csv');m=rd('figure_source_data/F2B_Blood_meta.csv').set_index('Program');l=rd('figure_source_data/F2C_Blood_leaveout.csv')
    for prog,d in b.groupby('Program'):
        close(meta(d['Hedges g'].to_numpy(),d.Variance.to_numpy()),m.loc[prog,['Summary Hedges g','95% CI low','95% CI high','Hartung–Knapp P']].to_numpy(float))
        for coh in d.Cohort:
            sub=d[d.Cohort.ne(coh)];close(meta(sub['Hedges g'].to_numpy(),sub.Variance.to_numpy())[0],l.loc[l.Program.eq(prog)&l['Omitted cohort'].eq(coh),'Summary Hedges g'].iloc[0])
    report['blood_meta_programmes']=7;report['blood_omission_estimates']=21
    m=rd('results/CSF_meta_11_proteins.csv')
    for _,r in m.iterrows():
        out=meta_two([r.discovery_beta,r.replication_beta],[r.discovery_se,r.replication_se])
        for key,val in out.items():close(val,r[key])
    for key,col in [('fixed_p','q_fixed_BH_11'),('random_p','q_random_BH_11'),('mhk_p','q_mhk_BH_11')]:close(bh(m[key]),m[col])
    assert ((m.mhk_low<0)&(m.mhk_high>0)).all();report['CSF_pooled_proteins']=11
    d=rd('followup_20260921/results/nerve_multiplicity_91.csv');close(bh(d.p_value),d.q_bh_91);assert len(d)==91
    d=rd('results/nerve_program_contrasts.csv');close(bh(d.p_mannwhitney),d.q_bh_9)
    a=rd('results/nerve_program_adjusted.csv');close(bh(a.p_hc3),a.q_bh_9_hc3)
    assert (d.q_bh_9>.05).all();report['nerve_expanded_tests']=91;report['new_nerve_contrasts']=9
    c=rd('source_data/nerve_clinical_plot_data.csv');expectedc=rd('results/nerve_clinical_correlations.csv').set_index('program')
    exact=[]
    for prog,d in c.groupby('program',sort=False):
        assert len(d)==9
        x=stats.rankdata(d.score);x-=x.mean();y=stats.rankdata(d.incat);y-=y.mean();norm=np.linalg.norm(x)*np.linalg.norm(y);rho=x@y/norm
        close(rho,expectedc.loc[prog,'rho']);extreme=0;total=0;it=itertools.permutations(y)
        while True:
            chunk=list(itertools.islice(it,20000))
            if not chunk:break
            scores=np.asarray(chunk)@x/norm;extreme+=int((np.abs(scores)>=abs(rho)-1e-12).sum());total+=len(chunk)
        assert total==math.factorial(9);p=extreme/total;close(p,expectedc.loc[prog,'exact_p']);exact.append((prog,p))
    old=rd('legacy_reference/analysis_update/source_data/CIDP_clinical.csv');q=bh(np.r_[old.p_value,[x[1] for x in exact]])[15:]
    for (prog,_),v in zip(exact,q):close(v,expectedc.loc[prog,'q_expanded_BH_17'])
    report['exact_clinical_tests']=2;report['permutations_per_test']=362880
    s=rd('figure_source_data/S2_Complement_count.csv').iloc[0];n=s.total_biopsies;x=s.positive_biopsies;p=x/n;z=stats.norm.ppf(.975);den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    close([p,mid-half,mid+half],[s.proportion,s.ci_low,s.ci_high]);report['published_count_check']='52/55; Wilson interval matches'
    f=rd('figure_source_data/S7_Normal_BNB_expression.csv');close(np.log2(f.fpkm+1),f.plotted_log2_fpkm_plus_1);assert len(f)==66
    status=json.loads((ROOT/'metadata/release.json').read_text());assert status['independent_nerve_validation'].startswith('not performed')
    report['independent_nerve_validation']='not performed; no disease mapping inferred'
    manifest=ROOT.parent/'MANIFEST_SHA256.json'
    if manifest.exists():
        hashes=json.loads(manifest.read_text())
        for name,h in hashes.items():assert hashlib.sha256((ROOT.parent/name).read_bytes()).hexdigest()==h,name
        report['manifest_files_matched']=len(hashes)
    print(json.dumps(report,indent=2));return report

if __name__=='__main__':main()
