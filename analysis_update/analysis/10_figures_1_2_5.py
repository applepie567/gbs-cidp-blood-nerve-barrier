from pathlib import Path
import io
import os
import json, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
for font_path in Path('/usr/local/share/fonts/truetype/timesnewroman').glob('*.TTF'):
    font_manager.fontManager.addfont(str(font_path))
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import ListedColormap, Normalize
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'results/figure_build'
SRC=ROOT/'source_data'
REPO=ROOT
OUT=ROOT
WORK.mkdir(parents=True,exist_ok=True)
FIG=OUT/'figures'
FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman','font.size':10,'text.color':'black',
 'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.titlesize':11,
 'axes.labelsize':10,'xtick.labelsize':9.5,'ytick.labelsize':9.5,'axes.linewidth':.7,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
TEAL='#158795';BLUE='#367AAD';ORANGE='#D98223';PURPLE='#8061A8';INK='#000000'
PANEL=15;TITLE=11
audit={}

def header(f,x,y,letter,title):
    f.text(x,y,letter,fontsize=PANEL,fontweight='bold',color=INK,va='baseline')
    f.text(x+.04,y,title,fontsize=TITLE,fontweight='bold',color=INK,va='baseline')

def clean(ax,grid=True):
    ax.spines[['top','right']].set_visible(False)
    if grid:ax.grid(axis='x',color='#dce1e4',lw=.55,zorder=0)
    ax.tick_params(length=3,pad=3)
    return ax

def save(f,n,name):
    f.canvas.draw()
    labels=[t for t in f.texts if t.get_text() in list('ABCDEF')]
    assert all(t.get_color()=='#000000' and t.get_fontsize()==PANEL for t in labels)
    audit[f'Figure_{n}']={'size_inches':list(f.get_size_inches()),'panel_font_pt':PANEL,
                        'panel_colour':'#000000','panel_labels':[t.get_text() for t in labels]}
    for ext in ['pdf','svg']:
        buf=io.BytesIO();f.savefig(buf,format=ext,dpi=400)
        (FIG/f'Figure_{n}_{name}.{ext}').write_bytes(buf.getvalue())
    import fitz
    pdf=fitz.open(FIG/f'Figure_{n}_{name}.pdf')
    for dpi,dest in [(400,FIG/f'Figure_{n}_{name}.png'),(160,WORK/f'preview_{n}.png')]:
        content=pdf[0].get_pixmap(dpi=dpi,alpha=False).tobytes('png')
        with open(dest,'wb') as fh:
            fh.write(content);fh.flush();os.fsync(fh.fileno())
        assert dest.stat().st_size==len(content)
    plt.close(f)

def figure1():
    f=plt.figure(figsize=(7.5,7.05));ax=f.add_axes([0,0,1,1]);ax.set_axis_off()
    header(f,.02,.967,'A','Data sources and analytical roles')
    cards=[(.035,.695,TEAL,'Acute GBS blood','Blood cohort comparison',
            'GSE211225 · whole blood\nGSE31014 · leukocytes\nPRJNA1293757 · PBMC'),
           (.52,.695,ORANGE,'GBS CSF','Published proteomic evidence',
            'PXD002911 · proteo-peptidomics\nComplement proteomics\nTMT and Olink studies'),
           (.035,.442,PURPLE,'CIDP peripheral nerve','Expression by donor and cell',
            'GSE285983 · sural nerve\nGSE107574 · published reference\nCIDP versus CIAP'),
           (.52,.442,BLUE,'CIDP genetics','Published genetic evidence',
            'GWAS in women\nMR and colocalization\nCandidate gene expression')]
    for x,y,c,title,sub,body in cards:
        w=.445;h=.225
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=.012',
                      linewidth=1.4,edgecolor=c,facecolor='white'))
        ax.add_patch(FancyBboxPatch((x+.004,y+h-.057),w-.008,.053,
                      boxstyle='round,pad=0.004,rounding_size=.008',linewidth=0,facecolor=c,alpha=.13))
        ax.text(x+.020,y+h-.029,title,va='center',fontsize=12,fontweight='bold',color=INK)
        ax.text(x+.020,y+.141,sub,va='center',fontsize=10.5,fontweight='bold',color=INK)
        ax.text(x+.020,y+.108,body,va='top',fontsize=10.5,linespacing=1.55,color=INK)
    header(f,.02,.380,'B','Complementary analyses')
    stages=[('Consistency','Recurrent\nblood effects','CXCL8 and\ncomplement\ngene modules',TEAL),
            ('CSF proteins','Independent\nCSF studies','Inflammation\nand structural\nproteolysis',ORANGE),
            ('Localization','Cell state and\ndonor level','Endothelial, immune\nand Schwann cell\nexpression',PURPLE),
            ('Candidate genes','Genetic\nevidence','Cellular expression\nof candidate\ngenes',BLUE)]
    for i,(title,sub,out,c) in enumerate(stages):
        x=.035+i*.244;w=.207
        ax.add_patch(FancyBboxPatch((x,.193),w,.151,boxstyle='round,pad=.006,rounding_size=.010',
                      linewidth=1.3,edgecolor=c,facecolor='white'))
        ax.text(x+w/2,.312,title,ha='center',va='center',fontsize=10.5,fontweight='bold')
        ax.text(x+w/2,.256,sub,ha='center',va='center',fontsize=10.5,linespacing=1.4)
        ax.add_patch(FancyBboxPatch((x,.025),w,.118,boxstyle='round,pad=.006,rounding_size=.010',
                      linewidth=1.1,edgecolor=c,facecolor=c,alpha=.13))
        ax.text(x+w/2,.084,out,ha='center',va='center',fontsize=10.5,linespacing=1.35)
        ax.add_patch(FancyArrowPatch((x+w/2,.182),(x+w/2,.154),arrowstyle='-|>',mutation_scale=10,lw=1.2,color=c))
        if i<3:ax.add_patch(FancyArrowPatch((x+w+.009,.268),(x+.244-.009,.268),
                  arrowstyle='-|>',mutation_scale=9,lw=1.0,color='#555555'))
    save(f,1,'study_architecture')

PROGRAMS=['CXCL8–CXCR1/2','LIF/OSM–gp130','Complement','Fc receptor',
          'Transendothelial migration','Interferon/JAK–STAT','Inflammatory monocyte']
PROGLABEL=['CXCL8–\nCXCR1/2','LIF/OSM–\ngp130','Complement','Fc receptor',
           'Transendothelial\nmigration','Interferon/\nJAK–STAT','Inflammatory\nmonocyte']

def figure2():
    d=pd.read_csv(SRC/'Blood_cohort_effects.csv');m=pd.read_csv(SRC/'Blood_meta.csv').set_index('Program').loc[PROGRAMS]
    cohorts=['Whole blood','Leukocytes','Monocytes']
    vals=d.pivot(index='Program',columns='Cohort',values='Hedges g').loc[PROGRAMS,cohorts].values
    fdr=d.pivot(index='Program',columns='Cohort',values='FDR').loc[PROGRAMS,cohorts].values
    f=plt.figure(figsize=(7.5,7.4))
    a=f.add_axes([.20,.20,.278,.745]);b=f.add_axes([.727,.640,.255,.305]);c=f.add_axes([.727,.20,.255,.305])
    header(f,.025,.965,'A','Cohort effects');header(f,.535,.965,'B','Cross-cohort synthesis')
    header(f,.535,.534,'C','Leave-one-cohort-out')
    # Keep black numerical labels legible while retaining the diverging scale.
    light_diverging=ListedColormap(plt.get_cmap('RdBu_r')(np.linspace(0,1,256))[:,:3]*.55+.45)
    im=a.imshow(vals,cmap=light_diverging,vmin=-4,vmax=4,aspect='auto')
    a.set_yticks(range(7),PROGLABEL);a.tick_params(axis='y',length=0,pad=5)
    a.set_xticks(range(3),['Whole\nblood','Leuko-\ncytes','Mono-\ncytes']);a.tick_params(axis='x',length=0,pad=7)
    for i in range(7):
        for j in range(3):a.text(j,i,f'{vals[i,j]:.2f}'+('*' if fdr[i,j]<.05 else ''),ha='center',va='center',
                                fontsize=10.5,fontweight='bold' if fdr[i,j]<.05 else 'normal',
                                color='black')
    cb=f.colorbar(im,cax=f.add_axes([.20,.072,.278,.019]),orientation='horizontal',ticks=[-4,-2,0,2,4]);cb.set_label('Hedges g',labelpad=2)
    y=np.arange(7);effect=m['Summary Hedges g'].values
    b.errorbar(effect,y,xerr=np.array([effect-m['95% CI low'].values,m['95% CI high'].values-effect]),fmt='o',
               color=TEAL,markersize=4.7,elinewidth=1.4,zorder=3)
    clean(b);b.axvline(0,color='#7a858a',lw=.8);b.set_xlim(-6.3,6.7);b.set_xticks([-6,-3,0,3,6]);b.set_ylim(6.5,-.5)
    b.set_yticks(y,PROGLABEL);b.set_xlabel('Hedges g (95% CI)',labelpad=4)
    loco=[]
    for prog in PROGRAMS:
        estimates=[]
        for omit in cohorts:
            g=d[(d.Program==prog)&(d.Cohort!=omit)];e=g['Hedges g'].values;v=g.Variance.values;w=1/v
            fixed=np.average(e,weights=w);q=np.sum(w*(e-fixed)**2)
            tau=max(0,(q-1)/(w.sum()-np.sum(w*w)/w.sum()))
            estimate=np.average(e,weights=1/(v+tau));estimates.append(estimate)
            loco.append({'Program':prog,'Omitted cohort':omit,'Summary Hedges g':estimate})
        i=PROGRAMS.index(prog);c.plot([min(estimates),max(estimates)],[i,i],lw=1,color='#a5adb2',zorder=1)
        for j,(e,marker,col) in enumerate(zip(estimates,['o','s','^'],[TEAL,BLUE,ORANGE])):
            c.scatter(e,i+[-.14,0,.14][j],marker=marker,color=col,s=24,edgecolor='white',linewidth=.4,zorder=4,
                      label=cohorts[j] if i==0 else None)
    clean(c);c.axvline(0,color='#7a858a',lw=.8);c.set_ylim(6.5,-.5);c.set_xlim(-1.3,1.72);c.set_xticks([-1,0,1])
    c.set_yticks(y,PROGLABEL);c.set_xlabel('Summary Hedges g',labelpad=4)
    handles,labels=c.get_legend_handles_labels()
    f.legend(handles,['Omit whole blood','Omit leukocytes','Omit monocytes'],loc='lower left',
             bbox_to_anchor=(.65,.017),frameon=False,fontsize=9.3,ncol=1,handletextpad=.35,labelspacing=.32)
    pd.DataFrame(loco).to_csv(SRC/'Figure_2C_leave_one_cohort_out.csv',index=False)
    audit['Figure_2_alignment']={'A_bottom':.20,'C_bottom':.20,'A_top':.945,'B_top':.945,
                                'B_height':.305,'C_height':.305}
    save(f,2,'acute_GBS_blood')



CELLS=['B_cell','Macrophage','BNB_EC','Pericyte_VSMC','Perineurium','Endoneurial_stroma','Epineurial_stroma','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
CELL_LABELS=['B cell','Macrophage','BNB EC','Pericyte/VSMC','Perineurium','Endoneurial stroma','Epineurial stroma','Myelinating SC','Nonmyelinating SC','Repair-damage SC']
GENES=['CDH4','DIRAS1','GNG7','SLC39A3']

def figure5():
    from redraw_figures_4_5 import make5
    make5()

if __name__ == "__main__":
    figure1(); figure2(); figure5()
