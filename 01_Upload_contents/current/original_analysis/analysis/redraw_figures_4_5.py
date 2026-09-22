"""Recreate the revised manuscript figures from the accompanying aggregate tables."""
from pathlib import Path
from io import BytesIO
import os
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'source_data'
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
for fp in Path('/usr/share').glob('fonts/**/NimbusRoman*.otf'):
    font_manager.fontManager.addfont(str(fp))
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Nimbus Roman', 'Times New Roman', 'DejaVu Serif'],
    'font.size': 9, 'text.color': 'black', 'axes.labelcolor': 'black',
    'xtick.color': 'black', 'ytick.color': 'black', 'axes.linewidth': .6,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white',
    'mathtext.fontset': 'stix',
})
BLUE, TEAL, ORANGE = '#367AAD', '#2A9D8F', '#E96B48'
GENES4 = ['CLDN5','ICAM1','VCAM1','LIFR','IL6ST','OSMR','CXCL8','C3','C3AR1','FCGR2A','JUN','SOX10','MPZ']
CELLS4 = ['Macrophage','BNB_EC','Pericyte','Perineurium','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
GENES5 = ['CDH4','DIRAS1','GNG7','SLC39A3']
CELLS5 = ['B_cell','Macrophage','BNB_EC','Pericyte_VSMC','Perineurium','Endoneurial_stroma','Epineurial_stroma','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
NAMES = {'B_cell':'B cell','Macrophage':'Macrophage','BNB_EC':'BNB EC',
    'Pericyte':'Pericyte',
    'Pericyte_VSMC':'Pericyte and VSMC','Perineurium':'Perineurium',
    'Endoneurial_stroma':'Endoneurial stroma','Epineurial_stroma':'Epineurial stroma',
    'Myelinating_SC':'Myelinating SC','Nonmyelinating_SC':'Nonmyelinating SC',
    'Repair_damage_SC':'Repair and damage SC','Other_EC':'Other EC',
    'T_NK':'T and NK cells','Granulocyte':'Granulocyte'}
SHORT4 = ['Macro','BNB EC','Pericyte','Perineurium','mySC','nmSC','repairSC']
SHORT5 = ['B cell','Macrophage','BNB EC','Pericyte and VSMC','Perineurium',
    'Endoneurial stroma','Epineurial stroma','mySC','nmSC','repairSC']
e4 = pd.read_csv(SRC/'Figure_4A_mean_expression.csv')
e5 = pd.read_csv(SRC/'Figure_5B_donor_mean_expression_and_detection.csv')

def save(fig, name):
    for ax in fig.axes:
        for collection in ax.collections:
            collection.set_rasterized(False)
    for ext in ['svg', 'pdf', 'png']:
        buf=BytesIO()
        fig.savefig(buf,format=ext,dpi=400)
        dest=OUT/f'{name}.{ext}'
        tmp=OUT/f'{name}.{ext}.tmp'
        with tmp.open('wb') as fh:
            fh.write(buf.getvalue());fh.flush();os.fsync(fh.fileno())
        tmp.replace(dest)
    plt.close(fig)

def title(fig, x, y, letter, label, fs=11):
    fig.text(x, y, letter, fontsize=15, weight='bold', va='baseline')
    fig.text(x+.039, y, label, fontsize=fs, weight='bold', va='baseline')

def plot4a(fig, rect, barrect, fs=8.4):
    ax=fig.add_axes(rect)
    mat=e4.pivot(index='gene',columns='cell_group',values='mean_log2cpm').loc[GENES4,CELLS4]
    xx,yy=np.meshgrid(np.arange(7),np.arange(13))
    sc=ax.scatter(xx.ravel(),yy.ravel(),c=mat.values.ravel(),s=27,cmap='viridis',vmin=-1,vmax=13,linewidths=0)
    ax.set(xlim=(-.5,6.5),ylim=(12.6,-.6))
    ax.set_xticks(range(7),SHORT4,rotation=43,ha='right',rotation_mode='anchor',fontsize=fs)
    ax.set_yticks(range(13),GENES4,fontsize=fs)
    for t in ax.get_yticklabels():t.set_fontstyle('italic')
    ax.tick_params(length=2,pad=3)
    ax.set_axisbelow(True)
    ax.grid(axis='y',color='#edf0f2',lw=.35)
    cb=fig.colorbar(sc,cax=fig.add_axes(barrect),orientation='horizontal',ticks=[-1,0,4,8,12])
    cb.ax.tick_params(labelsize=fs-.7,pad=2,length=2)
    cb.set_label(r'Mean donor log$_2$(CPM + 0.5)',fontsize=fs-.2,labelpad=2)
    return ax

def make4():
    fig=plt.figure(figsize=(7.5,6.4))
    title(fig,.025,.961,'A','Nerve cell expression')
    title(fig,.525,.961,'B','CIDP–CIAP module differences')
    plot4a(fig,[.129,.63,.313,.294],[.168,.514,.255,.011],8.1)
    # All retained panels are reconstructed from their unchanged plotted values.
    panels=['BNB_identity_integrity','Leukocyte_transmigration','CXCL8_CXCR1_2',
        'LIF_LIFR','Complement','Fc_receptor','Macrophage_state','Schwann_myelin_repair']
    cells=['Macrophage','BNB_EC','Perineurium','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
    data=pd.read_csv(SRC/'CIDP_module_effects.csv').query("comparison == 'CIDP_vs_CIAP'")
    v=data.pivot(index='panel',columns='cell_group',values='delta_mean_z').loc[panels,cells]
    q=data.pivot(index='panel',columns='cell_group',values='fdr_within_celltype_comparison').loc[panels,cells]
    ax=fig.add_axes([.68,.63,.241,.294])
    im=ax.pcolormesh(np.arange(7)-.5,np.arange(9)-.5,v,cmap='RdBu_r',vmin=-1.2,vmax=1.2,shading='flat',rasterized=False)
    ax.set(xlim=(-.5,5.5),ylim=(7.5,-.5))
    ax.set_yticks(range(8),['BNB identity','Migration','CXCL8','LIF and gp130','Complement','Fc receptor','Macrophage state','Schwann repair'],fontsize=8)
    ax.set_xticks(range(6),['Macro','BNB EC','Perineurium','mySC','nmSC','repairSC'],fontsize=8,rotation=43,ha='right',rotation_mode='anchor')
    ax.tick_params(length=0,pad=4)
    for i,j in zip(*np.where(q.values<.05)):
        ax.plot(j+.27,i-.27,marker='D',ms=2.8,color='#333333')
    cb=fig.colorbar(im,cax=fig.add_axes([.931,.63,.011,.294]),ticks=[-1,-.5,0,.5,1])
    cb.ax.tick_params(labelsize=7.3,length=2,pad=2)
    cb.set_label('Difference in module score',fontsize=7.6,labelpad=3)
    title(fig,.025,.429,'C','Cellular composition')
    title(fig,.525,.429,'D','Associations with INCAT')
    ax=fig.add_axes([.195,.083,.247,.312])
    fr=pd.read_csv(SRC/'CIDP_cell_fractions.csv').query("comparison == 'CIDP_vs_CIAP'").sort_values('difference_fraction',ascending=False)
    ax.scatter(fr.difference_fraction,range(len(fr)),s=18,color=BLUE,zorder=3)
    ax.set_yticks(range(len(fr)),[NAMES[x] for x in fr.cell_group],fontsize=7.9)
    ax.set(ylim=(len(fr)-.4,-.6),xlim=(-.035,.07))
    ax.set_xticks([-.02,0,.02,.04,.06])
    ax.set_xlabel('Difference in mean fraction\n(CIDP − CIAP)',fontsize=8.1,labelpad=3)
    d=fig.add_axes([.735,.083,.247,.312])
    cl=pd.read_csv(SRC/'CIDP_clinical.csv')
    cs={'BNB_EC':'EC','Macrophage':'Macro','Repair_damage_SC':'repairSC'}
    ps={'Complement':'Complement','Fc_receptor':'Fc receptor','LIF_LIFR':'LIF and gp130','Leukocyte_transmigration':'Transmigration','Schwann_myelin_repair':'Schwann repair'}
    for sig,col,lab in [(True,ORANGE,'P < 0.05'),(False,TEAL,'P ≥ 0.05')]:
        ix=np.where((cl.p_value<.05).to_numpy()==sig)[0]
        d.scatter(cl.iloc[ix].rho_INCAT,ix,s=18,color=col,label=lab,zorder=3)
    d.set_yticks(range(len(cl)),[cs[r.cell_group]+' · '+ps[r.panel] for r in cl.itertuples()],fontsize=7.3)
    d.set(ylim=(len(cl)-.35,-.65),xlim=(-1.02,1.05))
    d.set_xticks([-1,-.5,0,.5,1]);d.set_xlabel('Spearman ρ with INCAT',fontsize=8.1,labelpad=4)
    fig.legend(*d.get_legend_handles_labels(),loc='center',bbox_to_anchor=(.819,.412),ncol=2,frameon=False,fontsize=7,columnspacing=1,handlelength=.8,handletextpad=.3)
    for p in [ax,d]:
        p.spines[['top','right']].set_visible(False);p.axvline(0,color='#a2aab2',lw=.7,zorder=0)
        p.tick_params(axis='y',length=0,pad=3);p.tick_params(axis='x',length=2,pad=3,labelsize=8)
    save(fig,'Figure_4_CIDP_nerve')
    fig=plt.figure(figsize=(5.9,5.35))
    title(fig,.025,.948,'A','Nerve cell expression',12)
    plot4a(fig,[.20,.285,.735,.60],[.28,.093,.56,.018],10.4)
    save(fig,'Figure_4A_revised')

# Evidence labels are nominal categories, never a numerical or ordinal score.
EVIDENCE = [
 ['CXCL8\nrecruitment',
  ['Whole blood test','q = 0.045','Pooled P = 0.114'],
  ['Published findings','Innate immune and','inflammatory proteins'],
  ['Expression context','Myeloid genes and','endothelial trafficking']],
 ['Complement',
  ['Whole blood test','q = 0.0076','Pooled P = 0.234'],
  ['Published findings','Complement proteins'],
  ['Expression localisation','Macrophage complement genes']],
 ['Fc receptor',
  ['Whole blood test','q = 0.038','Mixed cohort directions'],
  ['Not assessed'],
  ['Macrophage test','q = 0.0087 within 8 modules','q = 0.087 across 80 tests']],
 ['Interferon\nresponse',
  ['Whole blood test','q = 0.045','Mixed cohort directions'],
  ['Published annotation','Antiviral defence'],
  ['Published atlas markers [27]','OAS1 and MX1 localisation','Published PNP versus control [27]','Lower MX1 in nmSC']],
 ['Endothelial\nadhesion',
  ['Migration module test','Whole blood q = 0.126','Mixed cohort directions'],
  ['Not assessed'],
  ['Expression localisation','ICAM1 and VCAM1 in BNB EC']],
 ['Schwann\nrepair',
  ['Not assessed'],
  ['Published findings','Axonal domain peptide','degradation'],
  ['Expression localisation','Schwann repair and stress genes']],
]

def evidence_map(fig,rect,fs=8.3):
    ax=fig.add_axes(rect);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    x=[0,.197,.454,.694,1.0]
    head=.155
    weights=np.array([1,1,1,1.40,1,1])
    heights=(1-head)*weights/weights.sum()
    # Header and all cells use the same neutral treatment for all evidence types.
    for c,lab in enumerate(['Programme','Blood','CSF','Peripheral nerve']):
        ax.add_patch(Rectangle((x[c],1-head),x[c+1]-x[c],head,facecolor='#e9edef',edgecolor='white',lw=1.3))
        ax.text((x[c]+x[c+1])/2,1-head/2,lab,ha='center',va='center',fontsize=fs+1.2,weight='bold')
    y=1-head
    for i,row in enumerate(EVIDENCE):
        h=heights[i];y-=h
        for c in range(4):
            ax.add_patch(Rectangle((x[c],y),x[c+1]-x[c],h,facecolor=('#f5f7f8' if i%2==0 else '#ffffff'),edgecolor='#dce1e4',lw=.45))
        ax.text(.012,y+h/2,row[0],ha='left',va='center',fontsize=fs+.4,linespacing=1.08)
        for c,lines in enumerate(row[1:]):
            n=len(lines)
            # Keep line spacing physical and symmetric within each cell.
            lineh=(h*.70)/(max(3,n)-1)
            start=y+h/2+lineh*(n-1)/2
            for j,txt in enumerate(lines):
                is_type=(j==0 or (i==3 and c==2 and j==2))
                txt=re.sub(r'\b(OAS1|MX1|ICAM1|VCAM1)\b',lambda m:r'$\mathit{'+m.group(0)+'}$',txt)
                ax.text((x[c+1]+x[c+2])/2,start-j*lineh,txt,ha='center',va='center',fontsize=fs if not is_type else fs-.05,weight='bold' if is_type else 'normal')
    return ax

def dot5b(fig,rect,fs=8.8):
    ax=fig.add_axes(rect)
    vals=e5.pivot(index='gene',columns='cell_group',values='mean_expression').loc[GENES5,CELLS5]
    pct=e5.pivot(index='gene',columns='cell_group',values='mean_percent_expressing').loc[GENES5,CELLS5]
    for i in range(4):
        ax.axhline(i,color='#e8ecef',lw=.45,zorder=0)
        for j in range(10):
            p=pct.iloc[i,j]
            if p==0:
                ax.plot(j,i,marker='x',ms=3,color='#7f898e',mew=.65,zorder=2)
            else:
                # Scatter area is directly proportional to the mean donor percentage.
                ax.scatter(j,i,s=6*p,c=[vals.iloc[i,j]],cmap='YlOrRd',vmin=0,vmax=1.6,edgecolors='#535a60',linewidths=.25,zorder=3)
    ax.set(xlim=(-.5,9.5),ylim=(3.55,-.55))
    ax.set_xticks(range(10),SHORT5,fontsize=fs-.5,rotation=39,ha='right',rotation_mode='anchor')
    ax.set_yticks(range(4),GENES5,fontsize=fs+.2)
    for t in ax.get_yticklabels():t.set_fontstyle('italic')
    for sp in ax.spines.values():sp.set_visible(False)
    ax.tick_params(length=0,pad=3)
    return ax

def dot5b_legends(fig,barrect,sizerect,fs=8):
    cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1.6),cmap='YlOrRd'),cax=fig.add_axes(barrect),orientation='horizontal',ticks=[0,.5,1,1.5])
    cb.ax.tick_params(labelsize=fs-.5,length=2,pad=2)
    cb.set_label(r'Mean donor log$_2$(CP10k + 1)',fontsize=fs,labelpad=2)
    ax=fig.add_axes(sizerect);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.5,.96,'Mean percentage of expressing nuclei',ha='center',va='top',fontsize=fs)
    xs=[.08,.26,.45,.66,.88]
    ax.plot(xs[0],.50,marker='x',ms=3,color='#7f898e',mew=.65)
    ax.text(xs[0],.04,'0%',ha='center',va='bottom',fontsize=fs-.4)
    for xx,p in zip(xs[1:],[1,5,20,50]):
        ax.scatter(xx,.50,s=6*p,color='#e9c984',edgecolors='#535a60',linewidths=.25)
        ax.text(xx,.04,f'{p}%',ha='center',va='bottom',fontsize=fs-.4)

def genetic_cards(fig,rect,fs=8.9):
    ax=fig.add_axes(rect);ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
    data=pd.read_csv(SRC/'CIDP_genetic_evidence.csv').set_index('gene')
    for i,gene in enumerate(GENES5):
        y=.765-i*.249
        ax.add_patch(FancyBboxPatch((0.003,y),.986,.213,boxstyle='round,pad=0.008,rounding_size=0.014',facecolor='#f5f7f8',edgecolor='#d6dde2',lw=.65))
        ax.text(.035,y+.112,gene,fontstyle='italic',weight='bold',fontsize=fs+1.2,va='center')
        row=data.loc[gene]
        if gene=='CDH4':
            line1='GWAS in women';line2='OR 4.37'
        else:
            line1='MR and colocalisation';line2=f'β {row.MR_beta:.2f}   PP.H4 {row.coloc_PP_H4:.2f}'
        ax.text(.385,y+.142,line1,weight='bold',fontsize=fs,va='center')
        ax.text(.385,y+.064,line2,fontsize=fs,va='center')

def difference5d(fig,rect,barrect,fs=8):
    ax=fig.add_axes(rect)
    data=pd.read_csv(SRC/'Genetic_CIDP_vs_CIAP.csv')
    # Column names are fixed by the archived source table.
    delta_col='delta_mean_expression'
    if delta_col not in data.columns:
        delta_col=next(c for c in data.columns if c.startswith('delta') or c.startswith('difference'))
    vals=data.pivot(index='cell_group',columns='gene',values=delta_col).loc[CELLS5,GENES5]
    im=ax.pcolormesh(np.arange(5)-.5,np.arange(11)-.5,vals,cmap='RdBu_r',vmin=-.8,vmax=.8,shading='flat',rasterized=False)
    ax.set(xlim=(-.5,3.5),ylim=(9.5,-.5))
    ax.set_yticks(range(10),[NAMES[x] for x in CELLS5],fontsize=fs)
    ax.set_xticks(range(4),GENES5,fontsize=fs,rotation=43,ha='right',rotation_mode='anchor')
    for t in ax.get_xticklabels():t.set_fontstyle('italic')
    ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,10,1),minor=True)
    ax.grid(which='minor',color='white',lw=.7)
    ax.tick_params(which='both',length=0,pad=3)
    for sp in ax.spines.values():sp.set_visible(False)
    cb=fig.colorbar(im,cax=fig.add_axes(barrect),orientation='horizontal',ticks=[-.8,0,.8])
    cb.ax.tick_params(labelsize=fs-.7,length=2,pad=2)
    return ax

def make5():
    fig=plt.figure(figsize=(7.5,9.8))
    title(fig,.023,.975,'A','Evidence types across compartments')
    evidence_map(fig,[.023,.622,.952,.332],7.65)
    title(fig,.023,.594,'B','CIDP nerve expression of genetic candidates')
    dot5b(fig,[.135,.465,.828,.107],8.6)
    dot5b_legends(fig,[.116,.337,.285,.009],[.482,.328,.482,.055],7.6)
    title(fig,.023,.280,'C','Published genetic evidence')
    title(fig,.535,.280,'D','CIDP–CIAP expression')
    genetic_cards(fig,[.023,.030,.442,.227],8.9)
    difference5d(fig,[.769,.080,.206,.177],[.769,.024,.206,.009],7.5)
    save(fig,'Figure_5_cross_compartment_genetics')
    fig=plt.figure(figsize=(7.5,4.20))
    title(fig,.023,.948,'A','Evidence types across compartments',12)
    evidence_map(fig,[.023,.038,.952,.855],8.6)
    save(fig,'Figure_5A_revised')
    fig=plt.figure(figsize=(7.5,4.4))
    title(fig,.023,.942,'B','CIDP nerve expression of genetic candidates',12)
    dot5b(fig,[.14,.455,.827,.385],10.2)
    dot5b_legends(fig,[.13,.133,.28,.020],[.49,.104,.478,.132],8.7)
    save(fig,'Figure_5B_revised')
    rows=[]
    for row in EVIDENCE:
        for comp,lines in zip(['Blood','CSF','Peripheral nerve'],row[1:]):
            rows.append({'programme':row[0].replace('\n',' '),'compartment':comp,
                'evidence_type':lines[0],'display_text':' | '.join(lines),
                'ordinal_evidence_rank':'not applicable'})
    pd.DataFrame(rows).to_csv(SRC/'Figure_5A_evidence_types.csv',index=False)

if __name__=='__main__':
    make4();make5()
    print('Saved revised Figures 4, 5, 4A, 5A and 5B in SVG, PDF and 400 dpi PNG.')
