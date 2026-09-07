from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
for font_path in Path('/usr/local/share/fonts/truetype/timesnewroman').glob('*.TTF'):
    font_manager.fontManager.addfont(str(font_path))

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / 'figures'
REPO = ROOT
plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 10.5,
    'text.color': 'black', 'axes.labelcolor': 'black', 'xtick.color': 'black',
    'ytick.color': 'black', 'axes.linewidth': .7, 'pdf.fonttype': 42,
    'savefig.facecolor': 'white'})
TEAL = '#158795'
BLUE = '#367AAD'
ORANGE = '#D98223'

def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='x', color='#dce1e4', lw=.55, zorder=0)
    ax.tick_params(length=3, pad=3)

def header(f,x,y,letter,title):
    f.text(x,y,letter,fontsize=15,fontweight='bold',color='#000000',va='baseline')
    f.text(x+.04,y,title,fontsize=11,fontweight='bold',color='#000000',va='baseline')

d = pd.read_csv(REPO/'results/tables/csf_kmezic_adjusted_meta.csv').set_index('protein')
d = d.loc[['SELE','IL2RA','THBD','CR2','CD38','CCL3','CD1C','NRP1','IL6','CD5']]
evidence = pd.read_csv(REPO/'results/tables/csf_published_evidence.csv')
programs = ['Protein accumulation','Structural proteolysis','Local inflammation',
            'Complement','Protein handling and neural structure']
coverage = evidence.groupby('program').agg(studies=('study','nunique'),platforms=('platform','nunique')).loc[programs]
coverage.to_csv(ROOT/'source_data/Figure_3C_coverage.csv')
f = plt.figure(figsize=(7.5,7.1))
a = f.add_axes([.11,.535,.35,.375])
b = f.add_axes([.615,.535,.30,.375])
c = f.add_axes([.43,.12,.50,.24])
header(f,.025,.952,'A','Adjusted cohort effects')
header(f,.535,.952,'B','Fixed effect synthesis')
header(f,.025,.388,'C','Coverage of complementary studies')
y = np.arange(len(d))
for name,offset,col,marker,label in [('discovery',.13,BLUE,'o','Discovery'),('replication',-.13,ORANGE,'s','Replication cohort')]:
    a.errorbar(d[name+'_beta'],y+offset,xerr=1.96*d[name+'_se'],fmt=marker,
               color=col,markersize=4.8,elinewidth=1.2,label=label,zorder=3)
clean(a)
a.axvline(0,color='#777777',lw=.8)
a.set_yticks(y,d.index,fontsize=11)
a.set_ylim(9.55,-.55)
a.set_xlim(-.8,4.1)
a.set_xticks([0,1,2,3,4])
a.set_xlabel('Adjusted coefficient (95% CI)',fontsize=10.5)
f.legend(*a.get_legend_handles_labels(),bbox_to_anchor=(.08,.431),loc='lower left',
         fontsize=9.8,ncol=2,frameon=False,columnspacing=.6,handlelength=1.0,handletextpad=.35)
b.errorbar(d.fixed_beta,y,xerr=np.array([d.fixed_beta-d.ci_low,d.ci_high-d.fixed_beta]),
           fmt='o',color=TEAL,markersize=5,elinewidth=1.45,zorder=3)
clean(b)
b.set_yticks(y,d.index,fontsize=11)
b.set_ylim(9.55,-.55)
b.set_xlim(-.35,2.65)
b.set_xticks([0,1,2])
b.set_xlabel('Pooled coefficient (95% CI)',fontsize=10.5)
b.text(1.12,1.025,'I²',transform=b.transAxes,ha='center',fontsize=11,fontweight='bold')
for i,v in enumerate(d.I2_percent):
    b.text(1.12,i,f'{v:.0f}%',transform=b.get_yaxis_transform(),ha='center',va='center',fontsize=10.5,color='black')
labels=['Protein accumulation','Structural proteolysis','Local inflammation',
        'Complement','Protein handling and\nneural structure']
for off,column,col,label in [(-.25,'studies',TEAL,'Studies'),(.25,'platforms',ORANGE,'Platforms')]:
    vals=coverage[column].to_numpy()
    c.barh(np.arange(5)+off,vals,height=.28,color=col,label=label,zorder=2)
    for i,v in enumerate(vals): c.text(v+.055,i+off,str(v),va='center',fontsize=10.5)
c.set_yticks(range(5),labels,fontsize=10.5)
c.set_ylim(4.7,-.7)
c.set_xlim(0,3.4)
c.set_xticks([0,1,2,3])
c.set_xlabel('Count')
clean(c)
f.legend(*c.get_legend_handles_labels(),bbox_to_anchor=(.52,.016),loc='lower left',
         fontsize=10,ncol=2,frameon=False,columnspacing=1,handlelength=1.3)
b.axvline(0,color='black',lw=.6)
f.savefig(WORK/'Figure_3_GBS_CSF.svg')
f.savefig(WORK/'Figure_3_GBS_CSF.pdf')
f.savefig(WORK/'Figure_3_GBS_CSF.png',dpi=400)

plt.close(f)
print('Figure 3 updated from existing coefficient and evidence tables')
