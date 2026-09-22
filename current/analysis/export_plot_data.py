"""Export only anonymous coordinates used in the two clinical scatter plots."""
import argparse
from pathlib import Path
import pandas as pd

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--private',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    s=pd.read_csv(a.private/'nerve_program_scores.csv')
    s=s.query('cell_group=="Repair_damage_SC" and disease=="CIDP"')[['program','score','incat']]
    assert len(s)==18
    s.sort_values(['program','score']).to_csv(a.out,index=False)
