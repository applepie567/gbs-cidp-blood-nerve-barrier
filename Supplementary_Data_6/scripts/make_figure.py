from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

R=Path(__file__).resolve().parents[1];D=R/'results';O=R/'outputs/Figures';O.mkdir(exist_ok=True,parents=True)
BLUE='#32688E';TEAL='#287E80';ORANGE='#C4652F';GREY='#6B7280'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
s=json.loads((D/'composition_summary.json').read_text())
fig,axes=plt.subplots(2,2,figsize=(8.5,6.6));fig.subplots_adjust(left=.105,right=.985,top=.93,bottom=.10,hspace=.62,wspace=.65)
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=26,fontweight='bold')
    ax.text(-.22,1.18,letter,transform=ax.transAxes,fontweight='bold',fontsize=13)
for j,scheme in enumerate(['donor','centre']):
    ax=axes[0,j];d=pd.read_csv(D/f'{scheme}_holdout_predictions.csv')
    panel(ax,'AB'[j],('Donor' if scheme=='donor' else 'Centre')+' held out')
    for centre,col in [('Münster',BLUE),('Essen',TEAL),('Würzburg',ORANGE)]:
        sub=d[d.centre==centre]
        ax.scatter(sub.nucleus_fraction,sub['index'],label=centre,c=col,s=25,edgecolor='white',lw=.4)
    rho=d['index'].corr(d.nucleus_fraction,method='spearman')
    ax.text(0,1.035,f'n = 37 donors, Spearman ρ = {rho:.3f}',transform=ax.transAxes,fontsize=8.4)
    ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    ax.set_xlabel('Observed macrophage nucleus fraction')
    ax.set_ylabel('Macrophage marker index')
    ax.legend(frameon=False,fontsize=7.8,loc='lower right',handletextpad=.4)
    ax.grid(alpha=.12)
ax=axes[1,0];d=pd.read_csv(D/'external_donor_scores.csv')
panel(ax,'C','Independent whole nerve cohort')
for diagnosis,col,label in [('CIDP',BLUE,'CIDP (n = 4)'),('VN',ORANGE,'Vasculitic neuropathy (n = 9)')]:
    sub=d[d.diagnosis==diagnosis]
    ax.scatter(sub.macrophage_index,sub.Fc_score,color=col,label=label,s=35,edgecolor='white',lw=.4)
ax.text(0,1.035,f'Spearman ρ = {s["external"]["rho_index_Fc"]:.3f}',transform=ax.transAxes,fontsize=8.4)
ax.set_xlabel('Macrophage marker index');ax.set_ylabel('Fc module score')
ax.legend(frameon=False,fontsize=7.3,loc='upper left',handletextpad=.35)
ax.grid(alpha=.12)
ax=axes[1,1];e=pd.read_csv(D/'external_conditional_models.csv')
panel(ax,'D','Fc difference after index adjustment')
ax.axvline(0,color='#999999',lw=.8,ls='--')
labels=['Unadjusted','Marker index\nadjusted','Mean probe\nsensitivity']
for j,row in e.iterrows():
    ax.errorbar(row.difference,j,xerr=[[row.difference-row.ci_low],[row.ci_high-row.difference]],fmt='o',color=[BLUE,TEAL,GREY][j],ms=5,capsize=3,lw=1.2)
ax.set_yticks(range(3),labels);ax.set_ylim(2.5,-.6)
ax.set_xlabel('CIDP − vasculitic neuropathy\nMean Fc score difference (95% CI)')
ax.grid(axis='x',alpha=.12)
for ext in ['png','pdf','svg','tiff']:
    fig.savefig(O/f'Supplementary_Figure_12.{ext}',dpi=350,bbox_inches='tight',pad_inches=.13)
plt.close(fig)
