"""Exact small-sample checks added after the initial exploratory correlations.

Enumerates all 9! assignments of INCAT ranks, retaining ties and their
multiplicities. No random resampling or statistical model selection is used.
The original 15 source P values remain the archived asymptotic Spearman tests.
The expanded BH family therefore explicitly mixes those existing tests with
the two exact new tests and is an exploratory multiplicity sensitivity.
"""
import argparse
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from stats_utils import bh


def run(private,out,legacy):
    s=pd.read_csv(private/'nerve_program_scores.csv').query('cell_group=="Repair_damage_SC" and disease=="CIDP"')
    rows=[]
    for program,d in s.groupby('program',sort=False):
        d=d.sort_values('sample')
        x=stats.rankdata(d.score);x=x-x.mean()
        y=stats.rankdata(d.incat);y=y-y.mean()
        norm=np.linalg.norm(x)*np.linalg.norm(y)
        observed=x@y/norm
        permutations=itertools.permutations(y)
        extreme=0;n=0
        while True:
            chunk=list(itertools.islice(permutations,20000))
            if not chunk:break
            values=np.array(chunk)@x/norm
            extreme+=int((np.abs(values)>=abs(observed)-1e-12).sum())
            n+=len(chunk)
        loo=[stats.spearmanr(d.drop(index=i).score,d.drop(index=i).incat).statistic for i in d.index]
        rows.append(dict(program=program,exact_p=extreme/n,permutations=n,
                         omission_rho_min=min(loo),omission_rho_max=max(loo)))
    current=pd.read_csv(out/'nerve_clinical_correlations.csv')
    current=current.drop(columns=[c for c in ['exact_p','permutations','omission_rho_min','omission_rho_max','q_exact_BH_2','q_expanded_BH_17'] if c in current])
    current=current.merge(pd.DataFrame(rows),on='program',validate='one_to_one')
    current['q_exact_BH_2']=bh(current.exact_p)
    old=pd.read_csv(legacy/'analysis_update/source_data/CIDP_clinical.csv')
    q=bh(np.r_[old.p_value.to_numpy(),current.exact_p.to_numpy()])
    current['q_expanded_BH_17']=q[15:]
    old['q_bh_with_two_exact_tests']=q[:15]
    current.to_csv(out/'nerve_clinical_correlations.csv',index=False)
    old.to_csv(out/'original_clinical_with_exact_extension.csv',index=False)
    print(current.to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['private','out','legacy']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();run(a.private,a.out,a.legacy)
