"""Rebuild strengthened figures from aggregate results and minimal plot data."""
import argparse
from io import BytesIO
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Bbox

BLUE='#367AAD'
ORANGE='#D96843'
TEAL='#298F87'
GREY='#687780'
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],
 'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'svg.fonttype':'none',
 'pdf.fonttype':42,'savefig.facecolor':'white','mathtext.fontset':'stix'})

def save(fig,out,name):
    for ext in ['png','pdf','svg']:
        buf=BytesIO()
        crop={'Figure_1':1.45,'Figure_3':.375,'Figure_5':.25}.get(name,0)
        bbox=Bbox.from_extents(0,crop,*fig.get_size_inches()) if crop else None
        fig.savefig(buf,format=ext,dpi=400,bbox_inches=bbox)
        tmp=out/f'{name}.{ext}.tmp';tmp.write_bytes(buf.getvalue());tmp.replace(out/f'{name}.{ext}')
    plt.close(fig)

def title(fig,x,y,letter,label):
    fig.text(x,y,letter,fontsize=14,fontweight='bold',va='baseline')
    fig.text(x+.038,y,label,fontsize=10,fontweight='bold',va='baseline')

def zero(ax):
    ax.axvline(0,color='#a7b1b7',lw=.7,zorder=0)
    ax.grid(axis='x',color='#edf0f2',lw=.5,zorder=0)
    ax.set_axisbelow(True)

def figure1(out):
    fig=plt.figure(figsize=(7.5,4.4))
    ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.5,.95,'Complementary immune profiles in inflammatory neuropathies',ha='center',fontsize=13,weight='bold')
    cards=[(.025,BLUE,'Blood','Acute GBS and healthy controls',
      'Whole blood and leukocytes\nMonocyte pseudobulk\n3 cohorts, 31 participants\n7 modules and cohort synthesis'),
      (.355,TEAL,'Cerebrospinal fluid','GBS, CIDP and healthy controls',
      'Age, sex and CCL19 adjusted\n11 paired GBS protein estimates\n16 disease comparison estimates\nSource tests and pooling uncertainty'),
      (.685,ORANGE,'Peripheral nerve','CIDP and axonal neuropathy',
      '9 CIDP and 11 CIAP donors\nMacrophages and endothelium\nSchwann cell programmes\nComposition and donor variation')]
    for x,c,h,sub,body in cards:
        ax.add_patch(FancyBboxPatch((x,.37),.29,.47,boxstyle='round,pad=0.01,rounding_size=.012',fc='#f6f8f9',ec='#d8e0e4',lw=.8))
        ax.add_patch(FancyBboxPatch((x,.74),.29,.10,boxstyle='round,pad=0.01,rounding_size=.012',fc=c,ec=c))
        ax.text(x+.145,.79,h,color='white',weight='bold',fontsize=12,ha='center',va='center')
        ax.text(x+.145,.68,sub,ha='center',weight='bold',fontsize=8.9)
        ax.text(x+.145,.515,body,ha='center',va='center',fontsize=8.2,linespacing=1.9)
    save(fig,out,'Figure_1')

def figure3(res,out):
    m=pd.read_csv(res/'CSF_meta_11_proteins.csv').set_index('protein')
    names=['IL8','SELE','IL2RA','CCL3','CR2','CD1C','THBD','NRP1','CD38','IL6','CD5']
    m=m.loc[names];y=np.arange(len(m))
    fig=plt.figure(figsize=(7.5,7.5))
    title(fig,.025,.958,'A','GBS versus healthy controls')
    title(fig,.52,.958,'B','Pooled GBS coefficients')
    a=fig.add_axes([.10,.585,.36,.32]);b=fig.add_axes([.59,.585,.36,.32])
    for prefix,col,off,label in [('discovery',BLUE,-.13,'Discovery'),('replication',ORANGE,.13,'Replication')]:
        a.errorbar(m[prefix+'_beta'],y+off,xerr=1.95996398454*m[prefix+'_se'],fmt='o',ms=3.3,color=col,lw=.8,label=label)
    for prefix,col,off,label in [('fixed',BLUE,-.13,'Fixed effect'),('random',ORANGE,.13,'Random effects')]:
        b.errorbar(m[prefix+'_beta'],y+off,xerr=[m[prefix+'_beta']-m[prefix+'_low'],m[prefix+'_high']-m[prefix+'_beta']],fmt='o',ms=3.3,color=col,lw=.8,label=label)
    for ax in [a,b]:
        zero(ax);ax.set_yticks(y,names);ax.set_ylim(10.6,-.7);ax.set_xlim(-1.3,4.2)
        ax.set_xlabel('Adjusted coefficient, log₂ NPX')
        ax.legend(frameon=False,loc='lower center',bbox_to_anchor=(.5,1.015),ncol=2,fontsize=8,handlelength=1.4,columnspacing=.8)
    title(fig,.025,.470,'C','GBS versus CIDP')
    title(fig,.52,.470,'D','CIDP versus healthy controls')
    d=pd.read_csv(res/'CSF_disease_comparisons.csv')
    for rect,contrast in [([.10,.115,.36,.30],'GBS versus CIDP'),([.59,.115,.36,.30],'CIDP versus HC')]:
        ax=fig.add_axes(rect);dd=d[d.contrast.eq(contrast)].reset_index(drop=True)
        for i,r in dd.iterrows():
            c=ORANGE if r.source_q<.05 else GREY
            ax.errorbar(r.beta,i,xerr=[[r.beta-r.normal_ci_low],[r.normal_ci_high-r.beta]],fmt='o',ms=4,color=c,lw=.9)
        zero(ax);ax.set_yticks(range(len(dd)),dd.protein);ax.set_ylim(len(dd)-.4,-.6)
        ax.set_xlabel('Adjusted coefficient, log₂ NPX');ax.set_xlim(-.4,3.7 if contrast.startswith('GBS') else 2.6)
    save(fig,out,'Figure_3')

def figure4(legacy,out):
    src=legacy/'analysis_update/source_data';fc=legacy/'strengthening_v41/results'
    fig=plt.figure(figsize=(7.5,6.9))
    title(fig,.025,.960,'A','Nerve cell expression')
    title(fig,.53,.960,'B','CIDP–CIAP module differences')
    genes=['CLDN5','ICAM1','VCAM1','LIFR','IL6ST','OSMR','CXCL8','C3','C3AR1','FCGR2A','JUN','SOX10','MPZ']
    cells=['Macrophage','BNB_EC','Pericyte_VSMC','Perineurium','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
    e=pd.read_csv(src/'Figure_4A_mean_expression.csv')
    # Cell labels in the original aggregate table are authoritative.
    if 'Pericyte_VSMC' not in e.cell_group.unique(): cells[2]='Pericyte'
    mat=e.pivot(index='gene',columns='cell_group',values='mean_log2cpm').loc[genes,cells]
    ax=fig.add_axes([.12,.60,.31,.31]);xx,yy=np.meshgrid(range(7),range(13))
    sc=ax.scatter(xx.ravel(),yy.ravel(),c=mat.values.ravel(),s=25,cmap='viridis',vmin=-1,vmax=13,lw=0)
    ax.set(xlim=(-.5,6.5),ylim=(12.6,-.6));ax.set_yticks(range(13),genes)
    for t in ax.get_yticklabels():t.set_fontstyle('italic')
    ax.set_xticks(range(7),['Macro','BNB EC','Pericyte','Perineurium','mySC','nmSC','repairSC'],rotation=43,ha='right',rotation_mode='anchor')
    ax.tick_params(length=2)
    cb=fig.colorbar(sc,cax=fig.add_axes([.16,.493,.26,.012]),orientation='horizontal',ticks=[0,4,8,12])
    cb.set_label('Mean donor log₂(CPM + 0.5)',fontsize=8,labelpad=2);cb.ax.tick_params(labelsize=7,length=2,pad=2)
    panels=['BNB_identity_integrity','Leukocyte_transmigration','CXCL8_CXCR1_2','LIF_LIFR','Complement','Fc_receptor','Macrophage_state','Schwann_myelin_repair']
    c2=['Macrophage','BNB_EC','Perineurium','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
    d=pd.read_csv(src/'CIDP_module_effects.csv').query("comparison == 'CIDP_vs_CIAP'")
    v=d.pivot(index='panel',columns='cell_group',values='delta_mean_z').loc[panels,c2]
    q=d.pivot(index='panel',columns='cell_group',values='fdr_within_celltype_comparison').loc[panels,c2]
    ax=fig.add_axes([.70,.60,.23,.31]);im=ax.imshow(v,aspect='auto',cmap='RdBu_r',vmin=-1.2,vmax=1.2)
    ax.set_yticks(range(8),['BNB identity','Migration','CXCL8','LIF and gp130','Complement','Fc receptor','Macrophage state','Myelin and repair'])
    ax.set_xticks(range(6),['Macro','BNB EC','Perineurium','mySC','nmSC','repairSC'],rotation=43,ha='right',rotation_mode='anchor');ax.tick_params(length=0)
    for i,j in zip(*np.where(q.values<.05)):ax.plot(j,i,'D',ms=3,color='#333333')
    cb=fig.colorbar(im,cax=fig.add_axes([.74,.493,.19,.012]),orientation='horizontal',ticks=[-1,0,1])
    cb.set_label('Difference in module score',fontsize=8,labelpad=2);cb.ax.tick_params(labelsize=7,length=2,pad=2)
    title(fig,.025,.393,'C','Cell fractions')
    title(fig,.53,.393,'D','Macrophage Fc receptor score')
    f=pd.read_csv(src/'CIDP_cell_fractions.csv')
    labels={'Macrophage':'Macrophage','T_NK':'T and NK','B_cell':'B cell','Bcell':'B cell','Granulocyte':'Granulocyte','BNB_EC':'BNB EC','Other_EC':'Other EC','Pericyte_VSMC':'Pericyte and VSMC','Perineurium':'Perineurium','Endoneurial_stroma':'Endoneurial stroma','Myelinating_SC':'mySC','Nonmyelinating_SC':'nmSC','Repair_damage_SC':'repairSC'}
    val='difference_fraction'
    ax=fig.add_axes([.22,.085,.23,.26]);ax.barh(range(len(f)),f[val],color=BLUE,height=.65)
    ax.set_yticks(range(len(f)),[labels.get(x,x) for x in f.cell_group],fontsize=7.3);ax.invert_yaxis();zero(ax)
    ax.set_xlabel('Mean fraction difference',fontsize=8)
    s=pd.read_csv(fc/'fc_sensitivity.csv')
    comp=pd.read_csv(fc/'fc_composition_adjustment.csv')
    r1=s[s.analysis.eq('All donors')].iloc[0] if 'All donors' in s.analysis.values else s.iloc[0]
    r2=s[s.analysis.eq('Age and center')].iloc[0]
    r3=comp[comp.model.str.contains('composition',case=False)].iloc[-1] if 'model' in comp else comp.iloc[-1]
    pd.DataFrame([r1,r2,r3]).to_csv(out/'Figure_4D_model_rows.csv',index=False)
    ax=fig.add_axes([.75,.125,.21,.22])
    for i,r in enumerate([r1,r2,r3]):
        lo=r.get('ci_low',r.get('ci95_low',np.nan));hi=r.get('ci_high',r.get('ci95_high',np.nan))
        ax.errorbar(r.delta,i,xerr=[[r.delta-lo],[hi-r.delta]],fmt='o',color=BLUE if i<2 else ORANGE,ms=4,lw=1)
    zero(ax);ax.set_yticks(range(3),['Unadjusted','Age and centre','Age, centre and\ncomposition'],fontsize=8)
    ax.set_ylim(2.65,-.65);ax.set_xlim(-1.25,.65);ax.set_xlabel('Score difference and 95% CI',fontsize=8)
    save(fig,out,'Figure_4')

def figure5(res,plotdata,out):
    d=pd.read_csv(res/'nerve_program_contrasts.csv');a=pd.read_csv(res/'nerve_program_adjusted.csv')
    fig=plt.figure(figsize=(7.5,6.6))
    title(fig,.025,.953,'A','Endothelial programmes')
    title(fig,.53,.953,'B','Schwann cell programmes')
    for rect,mask,lab in [([.22,.62,.23,.26],d.cell_group.eq('BNB_EC'),['Junctions','Transport','Adhesion']),
                          ([.77,.62,.19,.26],~d.cell_group.eq('BNB_EC'),['mySC myelin','mySC injury','nmSC myelin','nmSC injury','repairSC myelin','repairSC injury'])]:
        dd=d[mask].reset_index(drop=True);ax=fig.add_axes(rect)
        for i,r in dd.iterrows():
            ar=a[a.cell_group.eq(r.cell_group)&a.program.eq(r.program)].iloc[0]
            for rr,off,col in [(r,-.12,BLUE),(ar,.12,ORANGE)]:
                ax.errorbar(rr.delta,i+off,xerr=[[rr.delta-rr.ci_low],[rr.ci_high-rr.delta]],fmt='o',ms=3,color=col,lw=.8)
        zero(ax);ax.set_yticks(range(len(dd)),lab);ax.set_ylim(len(dd)-.4,-.6);ax.set_xlim(-1.9,1.8);ax.set_xticks([-1.5,0,1.5]);ax.set_xlabel('CIDP–CIAP score difference',fontsize=8)
    fig.text(.22,.528,'●  Unadjusted',color=BLUE,fontsize=9)
    fig.text(.48,.528,'●  Age and centre adjusted',color=ORANGE,fontsize=9)
    title(fig,.025,.425,'C','Myelin maintenance and disability')
    title(fig,.53,.425,'D','Injury response and disability')
    p=pd.read_csv(plotdata);c=pd.read_csv(res/'nerve_clinical_correlations.csv').set_index('program')
    for rect,prog,col in [([.115,.105,.325,.245],'Myelin maintenance',BLUE),([.615,.105,.325,.245],'Injury response',ORANGE)]:
        ax=fig.add_axes(rect);pp=p[p.program.eq(prog)];r=c.loc[prog]
        ax.scatter(pp.score,pp.incat,s=28,color=col,edgecolor='white',lw=.4)
        ax.set_xlabel(prog+' score',fontsize=9);ax.set_ylabel('INCAT disability',fontsize=9)
        ax.set_ylim(-.3,max(pp.incat)+.5);ax.set_yticks(sorted(pp.incat.unique()))
        ax.text(.02,1.045,f'ρ = {r.rho:.3f}   Exact P = {r.exact_p:.4f}   q₁₇ = {r.q_expanded_BH_17:.3f}',transform=ax.transAxes,fontsize=8.3)
    save(fig,out,'Figure_5')

def supplement5(legacy,out):
    src=legacy/'analysis_update/source_data'
    genes=['CDH4','DIRAS1','GNG7','SLC39A3']
    e=pd.read_csv(src/'Figure_5B_donor_mean_expression_and_detection.csv')
    cells=['B_cell','Macrophage','BNB_EC','Pericyte_VSMC','Perineurium','Endoneurial_stroma','Epineurial_stroma','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
    # Preserve the original source cell order if label conventions differ.
    if not set(cells).issubset(e.cell_group.unique()):cells=list(e.cell_group.drop_duplicates())
    labs={'Bcell':'B cell','B_cell':'B cell','Macrophage':'Macrophage','BNB_EC':'BNB EC','Pericyte_VSMC':'Pericyte and VSMC','Pericyte':'Pericyte and VSMC','Perineurium':'Perineurium','Endoneurial_stroma':'Endoneurial stroma','Epineurial_stroma':'Epineurial stroma','Myelinating_SC':'mySC','Nonmyelinating_SC':'nmSC','Repair_damage_SC':'repairSC'}
    v=e.pivot(index='gene',columns='cell_group',values='mean_expression').loc[genes,cells]
    pct=e.pivot(index='gene',columns='cell_group',values='mean_percent_expressing').loc[genes,cells]
    fig=plt.figure(figsize=(7.5,7.2));title(fig,.025,.96,'A','Expression of published genetic candidates')
    ax=fig.add_axes([.12,.68,.80,.22])
    for i,g in enumerate(genes):
        for j,c in enumerate(cells):
            p=pct.loc[g,c]
            if p==0:ax.plot(j,i,'x',ms=3,color=GREY)
            else:ax.scatter(j,i,s=5*p,c=[v.loc[g,c]],cmap='YlOrRd',vmin=0,vmax=1.6,lw=.25,edgecolor=GREY)
    ax.set_yticks(range(4),genes,fontstyle='italic');ax.set_xticks(range(10),[labs.get(c,c) for c in cells],rotation=35,ha='right',rotation_mode='anchor');ax.set_ylim(3.6,-.6);ax.set_xlim(-.5,9.5);ax.tick_params(length=0)
    cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1.6),cmap='YlOrRd'),cax=fig.add_axes([.12,.535,.30,.012]),orientation='horizontal');cb.set_label('Mean donor log₂(CP10k + 1)',fontsize=8);cb.ax.tick_params(labelsize=7,length=2,pad=2)
    leg=fig.add_axes([.54,.515,.40,.06]);leg.axis('off');leg.set(xlim=(0,1),ylim=(0,1))
    for x,p in zip([.10,.36,.65,.9],[1,5,20,50]):leg.scatter(x,.65,s=5*p,color='#e3a76d',edgecolor=GREY,lw=.25);leg.text(x,.02,f'{p}%',ha='center',fontsize=8)
    fig.text(.74,.58,'Mean expressing nuclei',ha='center',fontsize=8)
    title(fig,.025,.427,'B','Published association estimates')
    ax=fig.add_axes([.06,.11,.40,.27]);ax.axis('off')
    rows=[['CDH4','OR 4.37 (2.61–7.33)\nWomen, P = 1.49 × 10⁻⁸'],['DIRAS1','MR β −0.45, q = 0.00789\nPP.H4 0.94'],['GNG7','MR β 0.13, q = 0.00583\nPP.H4 0.90'],['SLC39A3','MR β −0.68, q = 0.00274\nPP.H4 0.90']]
    evidence=pd.read_csv(src/'CIDP_genetic_evidence.csv').set_index('gene')
    for row in rows[1:]:
        r=evidence.loc[row[0]];row[1]=f'MR β {r.MR_beta:.2f}, q = {r.MR_q_value:.4g}\nPP.H4 {r.coloc_PP_H4:.2f}'
    tab=ax.table(cellText=rows,colWidths=[.25,.75],cellLoc='left',bbox=[0,0,1,1]);tab.auto_set_font_size(False);tab.set_fontsize(8.7)
    for (i,j),cell in tab.get_celld().items():cell.set_edgecolor('#dbe2e6');cell.set_linewidth(.4);cell.set_facecolor('#f6f8f9' if i%2==0 else 'white')
    title(fig,.53,.427,'C','CIDP–CIAP expression differences')
    d=pd.read_csv(src/'Figure_5D_genetic_expression_contrasts.csv')
    column='delta_mean_expression' if 'delta_mean_expression' in d else 'delta'
    if column not in d:column='delta_mean'
    vv=d.pivot(index='gene',columns='cell_group',values=column).loc[genes,cells]
    ax=fig.add_axes([.62,.21,.32,.17]);im=ax.imshow(vv,cmap='RdBu_r',vmin=-.7,vmax=.7,aspect='auto')
    ax.set_yticks(range(4),genes,fontstyle='italic');ax.set_xticks(range(10),[labs.get(c,c) for c in cells],rotation=45,ha='right',fontsize=6.5,rotation_mode='anchor');ax.tick_params(length=0)
    cb=fig.colorbar(im,cax=fig.add_axes([.64,.065,.28,.01]),orientation='horizontal',ticks=[-.5,0,.5]);cb.set_label('Mean expression difference',fontsize=8);cb.ax.tick_params(labelsize=7,length=2,pad=2)
    save(fig,out,'Supplementary_Figure_5')

def supplement6(legacy,out):
    d=pd.read_csv(legacy/'analysis_update/source_data/CIDP_clinical.csv')
    fig=plt.figure(figsize=(7.5,5.0));ax=fig.add_axes([.42,.15,.51,.76])
    labels=[]
    cells={'BNB_EC':'BNB EC','Macrophage':'Macrophage','Repair_damage_SC':'Repair and damage SC'}
    panels={'LIF_LIFR':'LIF and gp130','Leukocyte_transmigration':'Migration','Complement':'Complement','Fc_receptor':'Fc receptor','Schwann_myelin_repair':'Myelin and repair'}
    for _,r in d.iterrows():labels.append(cells.get(r.cell_group,r.cell_group)+' · '+panels.get(r.panel,r.panel))
    ax.barh(range(len(d)),d.rho_INCAT,color=BLUE,height=.64);ax.set_yticks(range(len(d)),labels,fontsize=8);ax.invert_yaxis();zero(ax)
    ax.set_xlabel('Spearman correlation with INCAT disability');ax.set_xlim(-1,1)
    fig.text(.025,.97,'Original nerve module correlations',fontsize=12,weight='bold',va='top')
    fig.text(.42,.035,'All 15 correlations have q > 0.05.',fontsize=9)
    save(fig,out,'Supplementary_Figure_6')

def supplement7(res,out):
    d=pd.read_csv(res/'BNB_normal_reference_expression.csv')
    gene_col='gene';sample_col='sample';value_col='FPKM' if 'FPKM' in d else 'fpkm'
    mat=d.pivot(index=gene_col,columns=sample_col,values=value_col)
    genes=['CLDN5','OCLN','TJP1','CDH5','ABCB1','SLC1A1','MFSD2A','ICAM1','VCAM1','SELE','SELP']
    samples=list(d[sample_col].drop_duplicates());mat=mat.loc[genes,samples]
    fig=plt.figure(figsize=(6.5,5.4));ax=fig.add_axes([.20,.20,.63,.66]);im=ax.imshow(np.log2(mat+1),aspect='auto',cmap='viridis',vmin=0)
    ax.set_yticks(range(11),genes,fontstyle='italic');ax.set_xticks(range(6),['P3 culture','P8 culture','32P1','203P1','346P1','347P1'],rotation=35,ha='right',rotation_mode='anchor');ax.tick_params(length=0)
    ax.axvline(1.5,color='white',lw=2);ax.axhline(3.5,color='white',lw=1.5);ax.axhline(6.5,color='white',lw=1.5)
    cb=fig.colorbar(im,cax=fig.add_axes([.87,.20,.025,.66]));cb.set_label('log₂(FPKM + 1)')
    fig.text(.5,.945,'Normal human blood–nerve barrier reference',ha='center',weight='bold',fontsize=12)
    fig.text(.12,.035,'Two cultured endothelial preparations and four normal nerve preparations\nGSE107574 provides anatomical context without a disease comparison',fontsize=8.5,linespacing=1.4)
    save(fig,out,'Supplementary_Figure_7')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['legacy','results','plot-data','out']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    figure1(a.out);figure3(a.results,a.out);figure4(a.legacy,a.out)
    figure5(a.results,a.plot_data,a.out);supplement5(a.legacy,a.out)
    supplement6(a.legacy,a.out);supplement7(a.results,a.out)
