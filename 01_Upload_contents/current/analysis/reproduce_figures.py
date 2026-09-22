"""Regenerate all five main and seven supplementary figures from public tables."""
from pathlib import Path
import argparse, importlib.util, json, shutil, subprocess, sys, tempfile

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=Path('rebuilt_figures'))
    args=p.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    if out==ROOT/'figures':raise ValueError('Choose a separate output folder to preserve the manuscript figures.')
    py=sys.executable;script=ROOT/'analysis'
    subprocess.run([py,str(script/'make_figures.py'),'--legacy',str(ROOT/'legacy_reference'),'--results',str(ROOT/'results'),'--plot-data',str(ROOT/'source_data/nerve_clinical_plot_data.csv'),'--out',str(out)],check=True)
    # Separate processes preserve each original plot's rcParams.
    subprocess.run([py,str(script/'plot_figure2.py'),'--out',str(out)],check=True)
    subprocess.run([py,str(script/'plot_supplement1.py'),'--out',str(out)],check=True)
    subprocess.run([py,str(script/'plot_supplement2.py'),'--out',str(out)],check=True)
    subprocess.run([py,str(script/'plot_supplement34.py'),'--results',str(ROOT/'legacy_reference/strengthening_v41/results'),'--out',str(out)],check=True)
    names=[f'Figure_{i}' for i in range(1,6)]+[f'Supplementary_Figure_{i}' for i in range(1,8)]
    for name in names:
        for ext in ['png','pdf','svg']:
            assert (out/f'{name}.{ext}').stat().st_size>0,(name,ext)
    print(json.dumps({'figures':len(names),'formats':['png','pdf','svg'],'output':str(out),'note':'S2 is an editable reconstruction of the retained author-checked PNG; source observations are identical.'},indent=2))

if __name__=='__main__':main()
