from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
ROOT=Path(__file__).resolve().parents[1]
for p in Path('/usr/local/share/fonts/truetype/timesnewroman').glob('*.TTF'):font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Times New Roman','font.size':10,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','pdf.fonttype':42,'svg.fonttype':'none'})

def save(f,name):
 for t in f.findobj(Text):t.set_color('black');t.set_fontfamily('Times New Roman')
 for ext in ['png','pdf','svg']:f.savefig(ROOT/'figures'/f'{name}.{ext}',dpi=400,facecolor='white')
 plt.close(f)

d=pd.read_csv(ROOT/'results/tables/genetic_cell_localization_bootstrap.csv')
f,axes=plt.subplots(2,2,figsize=(9.2,9.4));f.subplots_adjust(left=.18,right=.97,top=.96,bottom=.09,hspace=.32,wspace=.65)
names={'BNB_EC':'BNB endothelium','B_cell':'B cell','Endoneurial_stroma':'Endoneurial stroma','Epineurial_stroma':'Epineurial stroma','Myelinating_SC':'Myelinating Schwann','Nonmyelinating_SC':'Nonmyelinating Schwann','Repair_damage_SC':'Repair/damage Schwann','Pericyte_VSMC':'Pericyte / VSMC','Other_EC':'Other endothelium','T_NK':'T / NK','Adipo':'Adipocyte','Mast_cell':'Mast cell'}
for ax,gene in zip(axes.flat,['CDH4','DIRAS1','GNG7','SLC39A3']):
 z=d[d.gene==gene].sort_values('mean_expression',ascending=False)
 ax.errorbar(z.mean_expression,range(len(z)),xerr=np.vstack([z.mean_expression-z.expression_ci_low,z.expression_ci_high-z.mean_expression]),fmt='o',color='#158795',ms=3.8,lw=1)
 ax.set_yticks(range(len(z)),[names.get(x,x) for x in z.cell_group],fontsize=8.5);ax.invert_yaxis();ax.set_title(gene,fontstyle='italic',fontweight='bold')
 ax.set_xlabel('Mean log2(CP10k + 1)');ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',color='#dddddd',lw=.5);ax.set_axisbelow(True)
save(f,'Supplementary_Figure_1_v39')

z=pd.read_csv(ROOT/'results/tables/cidp_external_nerve_validation.csv').iloc[0]
f=plt.figure(figsize=(7.5,4.5));ax=f.add_axes([.13,.53,.77,.22]);ax.errorbar(z['proportion']*100,0,xerr=[[100*(z['proportion']-z['ci_low'])],[100*(z['ci_high']-z['proportion'])]],fmt='o',color='#158795',ms=7,lw=2)
ax.set_xlim(75,100);ax.set_xticks([75,80,85,90,95,100]);ax.set_yticks([]);ax.spines[['top','left','right']].set_visible(False);ax.set_xlabel('Biopsies with endoneurial capillary C5b-9 deposition (%)')
f.text(.08,.90,'Published CIDP nerve evidence',fontsize=15,fontweight='bold')
f.text(.08,.80,'52/55 biopsies: 94.5% (Wilson 95% CI 85.1–98.1%)',fontsize=12)
f.text(.08,.34,'Additional source observations',fontsize=12,fontweight='bold')
f.text(.08,.27,'Proteomics: 9 CIDP cases and 2 controls',fontsize=10.5)
f.text(.08,.215,'Higher C3, C5 and C9; C6 elevated in only two CIDP cases',fontsize=10.5)
f.text(.08,.15,'CD68-positive macrophages and CD8-positive T cells',fontsize=11)
f.text(.08,.07,'Source: Stascheit et al. 2025; DOI: 10.1007/s00401-025-02936-w',fontsize=9.5)
save(f,'Supplementary_Figure_2_v39')
print('Supplementary figures regenerated from retained aggregate estimates')
