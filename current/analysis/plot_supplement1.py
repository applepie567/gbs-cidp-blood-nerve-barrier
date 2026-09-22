from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

def available_serif():
    for family in ['Times New Roman','Nimbus Roman','Liberation Serif','DejaVu Serif']:
        try:
            font_manager.findfont(font_manager.FontProperties(family=family),fallback_to_default=False)
            return family
        except ValueError:
            pass
    return 'serif'

from matplotlib.text import Text
ROOT=Path(__file__).resolve().parents[1]
import argparse
pa=argparse.ArgumentParser();pa.add_argument('--out',type=Path,required=True);OUT=pa.parse_args().out;OUT.mkdir(parents=True,exist_ok=True)
for p in Path('/usr/local/share/fonts/truetype/timesnewroman').glob('*.TTF'):font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':available_serif(),'font.size':10,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','pdf.fonttype':42,'svg.fonttype':'none'})

def save(f,name):
 for t in f.findobj(Text):t.set_color('black');t.set_fontfamily(available_serif())
 for ext in ['png','pdf','svg']:f.savefig(OUT/f'{name}.{ext}',dpi=400,facecolor='white')
 plt.close(f)

d=pd.read_csv(ROOT/'original_analysis/results/tables/genetic_cell_localization_bootstrap.csv')
f,axes=plt.subplots(2,2,figsize=(9.2,9.4));f.subplots_adjust(left=.18,right=.97,top=.96,bottom=.09,hspace=.32,wspace=.65)
names={'BNB_EC':'BNB endothelium','B_cell':'B cell','Endoneurial_stroma':'Endoneurial stroma','Epineurial_stroma':'Epineurial stroma','Myelinating_SC':'Myelinating Schwann','Nonmyelinating_SC':'Nonmyelinating Schwann','Repair_damage_SC':'Repair/damage Schwann','Pericyte_VSMC':'Pericyte / VSMC','Other_EC':'Other endothelium','T_NK':'T / NK','Adipo':'Adipocyte','Mast_cell':'Mast cell'}
for ax,gene in zip(axes.flat,['CDH4','DIRAS1','GNG7','SLC39A3']):
 z=d[d.gene==gene].sort_values('mean_expression',ascending=False)
 ax.errorbar(z.mean_expression,range(len(z)),xerr=np.vstack([z.mean_expression-z.expression_ci_low,z.expression_ci_high-z.mean_expression]),fmt='o',color='#158795',ms=3.8,lw=1)
 ax.set_yticks(range(len(z)),[names.get(x,x) for x in z.cell_group],fontsize=8.5);ax.invert_yaxis();ax.set_title(gene,fontstyle='italic',fontweight='bold')
 ax.set_xlabel('Mean log2(CP10k + 1)');ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',color='#dddddd',lw=.5);ax.set_axisbelow(True)
save(f,'Supplementary_Figure_1')

