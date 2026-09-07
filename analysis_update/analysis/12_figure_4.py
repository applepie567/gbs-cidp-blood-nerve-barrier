from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
for font_path in Path('/usr/local/share/fonts/truetype/timesnewroman').glob('*.TTF'):
    font_manager.fontManager.addfont(str(font_path))
from matplotlib.colors import ListedColormap
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / 'figures'
SRC = ROOT / 'source_data'
WORK.mkdir(exist_ok=True)

# Figure 4: retain the existing localization dots and reconstruct B-D from the
# exact plotted source tables. No biological estimates or significance tests change.
plt.rcParams.update({'font.family':'Times New Roman','font.size':8,
    'text.color':'black','axes.labelcolor':'black','xtick.color':'black',
    'ytick.color':'black','axes.linewidth':.6,'pdf.fonttype':42,
    'savefig.facecolor':'white'})
TEAL='#158795'; BLUE='#367AAD'; ORANGE='#E96B48'
f = plt.figure(figsize=(7.5, 5.575))

def title(x,y,letter,text):
    f.text(x,y,letter,fontsize=15,fontweight='bold',color='black',va='baseline')
    f.text(x+.04,y,text,fontsize=11,fontweight='bold',color='black',va='baseline')

title(.025,.950,'A','Nerve cell-state localization')
title(.530,.950,'B','CIDP–CIAP program shifts')
title(.025,.498,'C','Cellular composition')
title(.530,.498,'D','Associations with INCAT')
a = f.add_axes([.128,.641,.312,.286])
b = f.add_axes([.679,.641,.242,.286])
c = f.add_axes([.190,.088,.250,.374])
d = f.add_axes([.735,.088,.247,.374])

# Only the interior plot pixels are retained. All surrounding text is redrawn.
original=Image.open(ROOT/'analysis/assets/figure4_localization_reference.png').convert('RGB')
box=(249,130,812,572)
inside=original.crop(box)
a.imshow(inside,extent=(0,1,1,0),aspect='auto',interpolation='lanczos')
xcenters=np.array([288.0,357.4,426.8,496.2,565.6,635.0,704.4,773.8])
ycenters=150.5+np.arange(13)*(553.2-150.5)/12
a.set_xticks((xcenters-box[0])/(box[2]-box[0]),
    ['Macro','Gran','BNB EC','Pericyte','Perineurium','mySC','nmSC','repairSC'],
    fontsize=7.8,rotation=43,ha='right',rotation_mode='anchor')
a.set_yticks((ycenters-box[1])/(box[3]-box[1]),
    ['CLDN5','ICAM1','VCAM1','LIFR','IL6ST','OSMR','CXCL8','C3','C3AR1','FCGR2A','JUN','SOX10','MPZ'],fontsize=8)
a.tick_params(length=2,pad=3)

panels=['BNB_identity_integrity','Leukocyte_transmigration','CXCL8_CXCR1_2',
        'LIF_LIFR','Complement','Fc_receptor','Macrophage_state','Schwann_myelin_repair']
cells=['Macrophage','BNB_EC','Perineurium','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
effects=pd.read_csv(SRC/'CIDP_module_effects.csv')
effects=effects[effects.comparison=='CIDP_vs_CIAP']
values=effects.pivot(index='panel',columns='cell_group',values='delta_mean_z').loc[panels,cells]
q=effects.pivot(index='panel',columns='cell_group',values='fdr_within_celltype_comparison').loc[panels,cells]
original_palette=np.asarray(original)[133:571,1925,:][::-1]/255.
im=b.imshow(values,cmap=ListedColormap(original_palette),vmin=-1.2,vmax=1.2,aspect='auto')
b.set_yticks(range(8),['BNB identity','Migration','CXCL8','LIF/LIFR','Complement','Fc receptor','Macrophage state','Schwann repair'],fontsize=7.8)
b.set_xticks(range(6),['Macro','BNB EC','Perineurium','mySC','nmSC','repairSC'],fontsize=7.8,
             rotation=43,ha='right',rotation_mode='anchor')
b.tick_params(length=0,pad=4)
for i,j in zip(*np.where(q.values<.05)):
    b.plot(j+.27,i-.27,marker='D',markersize=2.4,color='#323232')
cb=f.colorbar(im,cax=f.add_axes([.931,.641,.010,.286]),ticks=[-1,-.5,0,.5,1])
cb.ax.tick_params(labelsize=7.2,length=2,pad=2)
cb.set_label('Δ module score',fontsize=7.5,labelpad=2)

fr=pd.read_csv(SRC/'CIDP_cell_fractions.csv')
fr=fr[fr.comparison=='CIDP_vs_CIAP'].sort_values('difference_fraction',ascending=False)
names={'Nonmyelinating_SC':'Nonmyelinating SC','Endoneurial_stroma':'Endoneurial stroma',
       'Perineurium':'Perineurium','Myelinating_SC':'Myelinating SC','BNB_EC':'BNB EC',
       'Macrophage':'Macrophage','B_cell':'B cell','T_NK':'T NK','Granulocyte':'Granulocyte',
       'Repair_damage_SC':'Repair damage SC','Other_EC':'Other EC','Pericyte_VSMC':'Pericyte'}
c.scatter(fr.difference_fraction,range(len(fr)),s=16,color=BLUE,zorder=3)
c.set_yticks(range(len(fr)),[names.get(x,x) for x in fr.cell_group],fontsize=7.8)
c.set_ylim(len(fr)-.4,-.6)
c.set_xlim(-.035,.070)
c.set_xticks([-.02,0,.02,.04,.06])
c.set_xlabel('Difference in mean fraction\n(CIDP − CIAP)',fontsize=8,labelpad=3,linespacing=1.2)

cl=pd.read_csv(SRC/'CIDP_clinical.csv')
cellshort={'BNB_EC':'EC','Macrophage':'Macro','Repair_damage_SC':'repairSC'}
pshort={'Complement':'Complement','Fc_receptor':'Fc receptor','LIF_LIFR':'LIF/LIFR',
        'Leukocyte_transmigration':'Transmigration','Schwann_myelin_repair':'Schwann repair'}
dl=[cellshort[r.cell_group]+' · '+pshort[r.panel] for r in cl.itertuples()]
for sig,col,lab in [(True,ORANGE,'Nominal P < 0.05'),(False,'#2A9D8F','P ≥ 0.05')]:
    ix=np.where((cl.p_value<.05).to_numpy()==sig)[0]
    d.scatter(cl.iloc[ix].rho_INCAT,ix,s=17,color=col,label=lab,zorder=3)
d.set_yticks(range(len(cl)),dl,fontsize=7.3)
d.set_ylim(len(cl)-.35,-.65)
d.set_xlim(-1.02,1.05)
d.set_xticks([-1,-.5,0,.5,1])
d.set_xlabel('Spearman ρ with INCAT',fontsize=8,labelpad=4)
f.legend(*d.get_legend_handles_labels(),loc='center',bbox_to_anchor=(.805,.479),
         ncol=2,frameon=False,fontsize=6.8,columnspacing=.8,
         handlelength=.75,handletextpad=.35,labelspacing=.35,borderpad=.1)
for ax in [c,d]:
    ax.spines[['top','right']].set_visible(False)
    ax.axvline(0,color='#a2aab2',lw=.7,zorder=0)
    ax.tick_params(axis='y',length=0,pad=3)
    ax.tick_params(axis='x',length=2,pad=3,labelsize=7.8)

f.savefig(WORK/'Figure_4_CIDP_nerve.png',dpi=400)
f.savefig(WORK/'Figure_4_CIDP_nerve.pdf')
f.savefig(WORK/'Figure_4_CIDP_nerve.svg')

plt.close(f)
print('Figure 4 axes and legend enlarged; source estimates retained')
