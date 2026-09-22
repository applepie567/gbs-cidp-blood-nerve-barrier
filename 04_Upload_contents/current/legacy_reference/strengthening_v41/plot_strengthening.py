from pathlib import Path
import argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

def draw(results,out):
    out.mkdir(parents=True,exist_ok=True)
    font_manager.fontManager.addfont('/usr/local/share/fonts/truetype/timesnewroman/Times.TTF')
    font_manager.fontManager.addfont('/usr/local/share/fonts/truetype/timesnewroman/Timesbd.TTF')
    plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    def forest(ax,d,labels,xlabel):
        y=np.arange(len(d))[::-1]
        ax.errorbar(d.delta,y,xerr=np.vstack([d.delta-d.ci_low,d.ci_high-d.delta]),fmt='o',color='#24566B',ecolor='#666666',capsize=3,markersize=5)
        ax.axvline(0,color='#999999',ls=':',lw=1);ax.set_yticks(y,labels);ax.set_xlabel(xlabel)
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0,pad=7)
        ax.set_ylim(-.8,len(d)-.2);ax.grid(axis='x',color='#EEEEEE',lw=.7);ax.set_axisbelow(True)
    def save(fig,name):
        for ext in ['pdf','svg','png','tiff']:
            extra={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
            fig.savefig(out/f'{name}.{ext}',dpi=400,bbox_inches='tight',facecolor='white',**extra)
        plt.close(fig)
    g=pd.read_csv(results/'fc_component_genes.csv')
    s=pd.read_csv(results/'fc_sensitivity.csv').iloc[:5].copy()
    c=pd.read_csv(results/'fc_composition_adjustment.csv');c['analysis']='Age, center and composition'
    s=pd.concat([s,c],ignore_index=True)
    slabels=['All donors','Male donors','Standardisation in 20 donors','Age and center','Age, center and sex','Age, center and composition']
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12.4,6.6),gridspec_kw={'width_ratios':[1,1.25]})
    forest(ax,g,[f'{r.gene}   q = {r.q_bh_available_genes:.3f}' for r in g.itertuples()],'Difference in mean log2(CPM + 0.5)')
    forest(bx,s,slabels,'Difference in mean module z scores')
    ax.set_title('A   Individual Fc module genes',loc='left',fontweight='bold',pad=16)
    bx.set_title('B   Donor sensitivity analyses',loc='left',fontweight='bold',pad=16)
    fig.subplots_adjust(left=.19,right=.99,bottom=.15,top=.90,wspace=.98)
    save(fig,'Supplementary_Figure_3_v41')
    a=pd.read_csv(results/'fc_annotation_sensitivity.csv')
    labels=[]
    for r in a.itertuples():
        q=r.q_eligible_refined_subtypes
        extra=f', q = {q:.3f}' if np.isfinite(q) else ''
        labels.append(f'{r.annotation} ({r.n_case} and {r.n_reference} donors{extra})')
    fig,ax=plt.subplots(figsize=(10.4,6.8))
    forest(ax,a,labels,'Difference in mean module z scores within each annotation')
    ax.set_title('Fc module scores under alternative macrophage annotations',loc='left',fontweight='bold',pad=18)
    ax.axhline(len(a)-3.5,color='#CCCCCC',lw=.8)
    fig.subplots_adjust(left=.50,right=.97,bottom=.14,top=.89)
    save(fig,'Supplementary_Figure_4_v41')
    print('Saved two figures in four formats')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();draw(a.results,a.out)
