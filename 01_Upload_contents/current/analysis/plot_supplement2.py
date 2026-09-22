"""Editable reconstruction of the author-checked S2; published inputs only."""
from pathlib import Path
import argparse
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    for fp in Path('/usr/local/share/fonts/truetype/timesnewroman').glob('*.TTF'):font_manager.fontManager.addfont(str(fp))
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],'font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
    z=pd.read_csv(ROOT/'figure_source_data/S2_Complement_count.csv').iloc[0]
    f=plt.figure(figsize=(7.5,4.5));ax=f.add_axes([.12,.56,.80,.18])
    ax.errorbar(100*z.proportion,0,xerr=[[100*(z.proportion-z.ci_low)],[100*(z.ci_high-z.proportion)]],fmt='o',color='#158795',ms=6,lw=2)
    ax.set_xlim(75,100);ax.set_xticks([75,80,85,90,95,100]);ax.set_yticks([]);ax.spines[['top','left','right']].set_visible(False)
    ax.set_xlabel('Biopsies with endoneurial capillary C5b-9 deposition (%)',fontsize=10,labelpad=7)
    f.text(.07,.94,'Published CIDP nerve observations',fontsize=17,fontweight='bold')
    f.text(.07,.855,f'{int(z.positive_biopsies)} of {int(z.total_biopsies)} biopsies  {100*z.proportion:.1f}% (Wilson 95% CI {100*z.ci_low:.1f}–{100*z.ci_high:.1f}%)',fontsize=11.5)
    f.text(.07,.40,'Additional source observations',fontsize=13,fontweight='bold')
    texts=['Proteomic subset:  9 CIDP cases and 2 controls','Higher C3, C5 and C9 abundance was reported.','C6 protein was elevated in 2 of 9 cases.','Histology identified CD68 positive macrophages and CD8 positive T cells.']
    for y,text in zip([.32,.25,.18,.11],texts):f.text(.07,y,text,fontsize=10.5)
    f.text(.07,.035,'Stascheit et al. 2025   https://doi.org/10.1007/s00401-025-02936-w',fontsize=9)
    for ext in ['png','pdf','svg']:f.savefig(a.out/f'Supplementary_Figure_2.{ext}',dpi=400,facecolor='white')
    plt.close(f)

if __name__=='__main__':main()
