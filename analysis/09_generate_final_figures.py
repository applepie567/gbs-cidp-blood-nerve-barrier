from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from scipy import stats

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "recovered/package/Additional_file_1_machine_readable_results.xlsx"
OUT = ROOT / "output/figures"
RAW = ROOT / "output/figure_source_data"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

NAVY = "#183B56"
BLUE = "#2E6F9E"
TEAL = "#2A9D8F"
ORANGE = "#E76F51"
GOLD = "#E9C46A"
PURPLE = "#7A6FA8"
GREY = "#6B7280"
LIGHT = "#E8EEF3"
INK = "#20262E"
PALETTE = [BLUE, TEAL, ORANGE, GOLD, PURPLE, NAVY]

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "axes.edgecolor": "#4B5563",
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 600,
})

MODULE_MAP = {
    "CXCL8_CXCR1_2": "CXCL8–CXCR1/2",
    "CXCL8_axis": "CXCL8–CXCR1/2",
    "LIF_LIFR_IL6ST": "LIF/OSM–gp130",
    "IL6_family_LIF_axis": "LIF/OSM–gp130",
    "Complement": "Complement",
    "Fc_receptor": "Fc receptor",
    "Transendothelial_migration": "Transendothelial migration",
    "Interferon_JAK_STAT": "Interferon/JAK–STAT",
    "Inflammatory_monocyte": "Inflammatory monocyte",
}
MODULE_ORDER = ["CXCL8–CXCR1/2", "Inflammatory monocyte", "Complement", "Transendothelial migration", "Interferon/JAK–STAT", "LIF/OSM–gp130", "Fc receptor"]

def read(sheet):
    return pd.read_excel(SOURCE, sheet_name=sheet)

def panel(ax, letter, title):
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", color=NAVY, va="bottom")
    ax.set_title(title, loc="left", pad=10, color=NAVY)

def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=.7, color="#4B5563")

def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

def write(df, name):
    df.to_csv(RAW / f"{name}.csv", index=False)

def variance_g(g, n1, n0):
    return (n1+n0)/(n1*n0) + g*g/(2*max(n1+n0-2, 1))

def random_meta(sub):
    y=sub.effect.to_numpy(float); v=sub.variance.to_numpy(float); k=len(y)
    wf=1/v; fixed=(wf*y).sum()/wf.sum(); q=(wf*(y-fixed)**2).sum()
    c=wf.sum()-(wf**2).sum()/wf.sum(); tau=max(0,(q-(k-1))/c) if c>0 else 0
    w=1/(v+tau); mu=(w*y).sum()/w.sum(); se=math.sqrt(1/w.sum())
    crit=stats.t.ppf(.975,max(k-1,1)); lo=mu-crit*se; hi=mu+crit*se
    i2=max(0,(q-(k-1))/q*100) if q>0 else 0
    p=2*stats.t.sf(abs(mu/se),max(k-1,1)) if se else 1
    return pd.Series(dict(effect=mu,low=lo,high=hi,p=p,tau2=tau,I2=i2,k=k))

def collect_gbs_effects():
    rows=[]
    d=read("GSE211225 Modules")
    for _,r in d[d.comparison.eq("Acute_GBS_vs_HC")].iterrows(): rows.append(["Whole blood",MODULE_MAP[r.module],r.hedges_g,r.n_case,r.n_reference,r.module_fdr_within_comparison])
    d=read("GSE31014 Modules")
    for _,r in d.iterrows(): rows.append(["Leukocytes",MODULE_MAP[r.module],r.hedges_g,r.n_case,r.n_reference,r.module_fdr_within_comparison])
    d=read("GBS Module Effects")
    for _,r in d[d.cell_type.eq("CD11b")].iterrows(): rows.append(["CD11b+",MODULE_MAP[r.module],r.hedges_g,r.n_gbs,r.n_hc,r.module_fdr_within_cell_type])
    d=read("PRJNA Modules")
    for _,r in d[d.compartment.eq("Marker_defined_monocyte")].iterrows(): rows.append(["Monocytes",MODULE_MAP[r.module],r.hedges_g,r.n_case,r.n_reference,r.module_fdr_within_compartment])
    out=pd.DataFrame(rows,columns=["cohort","module","effect","n_case","n_control","fdr"])
    out["variance"]=[variance_g(*x) for x in out[["effect","n_case","n_control"]].itertuples(index=False,name=None)]
    return out

def figure1():
    design=read("Dataset Design")
    write(design,"Figure_1A_dataset_design")
    fig,ax=plt.subplots(figsize=(11.2,5.8)); ax.set_xlim(0,12); ax.set_ylim(0,6); ax.axis("off")
    label_x=.28; title_x=.70
    ax.text(label_x,5.65,"A",fontsize=14,fontweight="bold",color=NAVY)
    ax.text(title_x,5.65,"Study architecture",fontsize=13,fontweight="bold",color=NAVY)
    centers=[1.55,4.50,7.45,10.40]
    stages=[("Acute GBS blood","4 cohorts\nwhole blood · leukocytes\nCD11b+ · PBMC",BLUE),
            ("Longitudinal serum","20 paired patients\nacute → 1 year",TEAL),
            ("CIDP target nerve","37 donors · 365,708 nuclei\nmacrophage · BNB\nSchwann-cell states",ORANGE),
            ("Mechanistic synthesis","recruitment → barrier\neffector → repair",PURPLE)]
    box_w=2.28; box_y=3.18; box_h=1.58
    for (title,txt,color),x in zip(stages,centers):
        box=FancyBboxPatch((x-box_w/2,box_y),box_w,box_h,boxstyle="round,pad=.06,rounding_size=.12",fc="white",ec=color,lw=1.8)
        ax.add_patch(box); ax.text(x,4.36,title,ha="center",va="center",fontweight="bold",color=color,fontsize=9.35 if title=="Mechanistic synthesis" else 10)
        ax.text(x,3.72,txt,ha="center",va="center",color=INK,fontsize=8.5,linespacing=1.35)
    for x1,x2 in zip(centers[:-1],centers[1:]):
        ax.add_patch(FancyArrowPatch((x1+box_w/2+.10,3.96),(x2-box_w/2-.10,3.96),arrowstyle="-|>",mutation_scale=13,lw=1.3,color=GREY))
    programs=[("CXCL8 recruitment",BLUE),("Complement/Fc",TEAL),("BNB adhesion",GOLD),("Macrophage states",ORANGE),("Schwann repair",PURPLE)]
    ax.text(label_x,2.45,"B",fontsize=14,fontweight="bold",color=NAVY); ax.text(title_x,2.45,"Prespecified biological programs",fontsize=13,fontweight="bold",color=NAVY)
    for i,(t,c) in enumerate(programs):
        x=.7+i*2.27; ax.add_patch(FancyBboxPatch((x,.85),1.95,.72,boxstyle="round,pad=.04,rounding_size=.09",fc=c+"18",ec=c,lw=1.2)); ax.text(x+.975,1.21,t,ha="center",va="center",fontsize=8.5,color=INK)
    save(fig,"Figure_1_study_architecture")

def figure2():
    effects=collect_gbs_effects(); write(effects,"Figure_2A_cohort_module_effects")
    meta=effects.groupby("module",sort=False).apply(random_meta,include_groups=False).reset_index(); meta["module"]=pd.Categorical(meta.module,MODULE_ORDER,ordered=True); meta=meta.sort_values("module"); write(meta,"Figure_2B_random_effects_meta")
    fig=plt.figure(figsize=(12.5,8.2)); gs=fig.add_gridspec(2,2,width_ratios=[1.05,1],height_ratios=[1,1],wspace=.62,hspace=.48)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","Cross-cohort acute GBS effects")
    piv=effects.pivot(index="module",columns="cohort",values="effect").reindex(MODULE_ORDER)[["Whole blood","Leukocytes","CD11b+","Monocytes"]]
    im=ax.imshow(piv,cmap=mpl.colors.LinearSegmentedColormap.from_list("div",["#2C7BB6","#F7F7F7","#D7191C"]),vmin=-3,vmax=3,aspect="auto")
    ax.set_xticks(range(4),piv.columns,rotation=25,ha="right"); ax.set_yticks(range(len(piv)),piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v=piv.iloc[i,j]; ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=7.6,color="white" if abs(v)>1.8 else INK)
    cb=fig.colorbar(im,ax=ax,fraction=.036,pad=.025); cb.set_label("Hedges g",fontsize=8)
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Random-effects synthesis")
    y=np.arange(len(meta))[::-1]
    for yi,(_,r) in zip(y,meta.iterrows()):
        ax.plot([r.low,r.high],[yi,yi],color=GREY,lw=1.5,zorder=1); ax.scatter(r.effect,yi,s=42,color=TEAL if r.p<.05 else NAVY,edgecolor="white",lw=.7,zorder=2)
    ax.axvline(0,color="#9CA3AF",lw=.9); ax.set_yticks(y,meta.module); ax.set_xlabel("Random-effects Hedges g (95% CI)"); clean(ax)
    ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=TEAL,markeredgecolor="white",label="Nominal P < 0.05"),
                       Line2D([0],[0],marker="o",color="none",markerfacecolor=NAVY,markeredgecolor="white",label="P ≥ 0.05")],
              frameon=False,loc="upper left",fontsize=7.2,borderaxespad=.2)
    d=read("GSE211225 Modules"); d=d[d.comparison.isin(["Acute_GBS_vs_HC","Acute_GBS_vs_Postacute_GBS"])].copy(); d["module_label"]=d.module.map(MODULE_MAP); write(d,"Figure_2C_whole_blood_phase_effects")
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","Whole-blood phase contrasts")
    for k,(comp,col,off) in enumerate([("Acute_GBS_vs_HC",BLUE,-.13),("Acute_GBS_vs_Postacute_GBS",ORANGE,.13)]):
        s=d[d.comparison.eq(comp)].set_index("module_label").reindex(MODULE_ORDER); yy=np.arange(len(s))[::-1]+off
        ax.scatter(s.hedges_g,yy,s=34,color=col,label="Acute vs HC" if k==0 else "Acute vs post-acute",edgecolor="white",lw=.5)
    ax.set_xlim(left=0); ax.set_yticks(np.arange(7)[::-1],MODULE_ORDER); ax.set_xlabel("Hedges g"); ax.legend(frameon=False,loc="upper left",fontsize=7.5,borderaxespad=.2); clean(ax)
    d=read("PRJNA Sample Scores"); d=d[d.compartment.eq("Marker_defined_monocyte")].copy(); d["module_label"]=d.module.map(MODULE_MAP); write(d,"Figure_2D_monocyte_sample_scores")
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","Monocyte localization in PBMCs")
    sub=d[d.module_label.isin(["CXCL8–CXCR1/2","Inflammatory monocyte","Fc receptor"])]
    xpos={m:i for i,m in enumerate(["CXCL8–CXCR1/2","Inflammatory monocyte","Fc receptor"])}
    rng=np.random.default_rng(11)
    for cond,c,shift in [("HC",BLUE,-.12),("Acute_GBS",ORANGE,.12)]:
        for m,g in sub[sub.condition.eq(cond)].groupby("module_label"):
            x=xpos[m]+shift+rng.normal(0,.025,len(g)); ax.scatter(x,g.score_z,s=30,color=c,alpha=.9,edgecolor="white",lw=.5)
            ax.plot([xpos[m]+shift-.07,xpos[m]+shift+.07],[g.score_z.mean()]*2,color=c,lw=2)
    ax.axhline(0,color="#9CA3AF",lw=.8); ax.set_xticks(range(3),["CXCL8", "Inflammatory\nmonocyte","Fc receptor"]); ax.set_ylabel("Sample module score (z)"); ax.set_ylim(-.5,.68); clean(ax)
    ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=BLUE,label="HC"),Line2D([0],[0],marker="o",color="none",markerfacecolor=ORANGE,label="Acute GBS")],frameon=False,fontsize=7.5,loc="upper left",borderaxespad=.2)
    fig.suptitle("Acute GBS converges on a myeloid recruitment program",x=.02,y=.995,ha="left",fontsize=15,fontweight="bold",color=NAVY)
    save(fig,"Figure_2_acute_GBS_crosscohort")

def figure3():
    scores=read("GBS Proteome Sample Scores"); scores["module_label"]=scores.module.map(MODULE_MAP); write(scores,"Figure_3A_proteome_sample_scores")
    eff=read("GBS Proteome Modules"); eff["module_label"]=eff.module.map(MODULE_MAP); write(eff,"Figure_3B_proteome_module_effects")
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,2,wspace=.42,hspace=.47)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","Paired acute-to-one-year trajectories")
    mods=["CXCL8–CXCR1/2","Interferon/JAK–STAT","LIF/OSM–gp130"]
    paired=scores[scores.module_label.isin(mods)&scores.condition.isin(["GBS_Acute","GBS_Recovery"])].pivot_table(index=["matching_ID","module_label"],columns="condition",values="score_z").dropna().reset_index()
    write(paired,"Figure_3A_paired_scores")
    for i,m in enumerate(mods):
        s=paired[paired.module_label.eq(m)]
        for _,r in s.iterrows(): ax.plot([i-.16,i+.16],[r.GBS_Acute,r.GBS_Recovery],color="#C7D2DA",lw=.7,alpha=.8)
        ax.scatter(np.full(len(s),i-.16),s.GBS_Acute,s=18,color=ORANGE,zorder=2); ax.scatter(np.full(len(s),i+.16),s.GBS_Recovery,s=18,color=TEAL,zorder=2)
    ax.axhline(0,color="#9CA3AF",lw=.8); ax.set_xticks(range(3),["CXCL8","Interferon","LIF/OSM"]); ax.set_ylabel("Protein-program score (z)"); clean(ax)
    ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=ORANGE,label="Acute"),Line2D([0],[0],marker="o",color="none",markerfacecolor=TEAL,label="One year")],frameon=False,fontsize=7.5)
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","Within-patient recovery effects")
    s=eff[eff.comparison.eq("Acute_GBS_vs_Recovery_1y")].set_index("module_label").reindex(MODULE_ORDER).reset_index(); y=np.arange(len(s))[::-1]
    ax.axvline(0,color="#9CA3AF",lw=.8); ax.scatter(s.dz,y,c=[TEAL if q<.05 else NAVY for q in s.paired_t_p],s=42,edgecolor="white",lw=.6)
    for yi,(_,r) in zip(y,s.iterrows()): ax.text(r.dz+.06,yi,f"P={r.paired_t_p:.3f}",va="center",fontsize=7,color=GREY)
    ax.set_yticks(y,s.module_label); ax.set_xlabel("Paired standardized effect, dz"); clean(ax)
    ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=TEAL,markeredgecolor="white",label="Paired P < 0.05"),
                       Line2D([0],[0],marker="o",color="none",markerfacecolor=NAVY,markeredgecolor="white",label="Paired P ≥ 0.05")],
              frameon=False,fontsize=7.2,loc="upper left",borderaxespad=.2)
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","Cross-sectional protein effects")
    comps=["Acute_GBS_vs_HC","Recovery_1y_vs_HC"]; cols=[ORANGE,TEAL]
    for k,(comp,c) in enumerate(zip(comps,cols)):
        z=eff[eff.comparison.eq(comp)].set_index("module_label").reindex(MODULE_ORDER); yy=np.arange(7)[::-1]+(-.12 if k==0 else .12); ax.scatter(z.hedges_g,yy,s=34,color=c,edgecolor="white",lw=.5,label="Acute vs HC" if k==0 else "One year vs HC")
    ax.axvline(0,color="#9CA3AF",lw=.8); ax.set_yticks(np.arange(7)[::-1],MODULE_ORDER); ax.set_xlabel("Hedges g"); ax.legend(frameon=False,fontsize=7.5,loc="lower right"); clean(ax)
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","Phase-specific interpretation")
    ax.axis("off")
    boxes=[("Recruitment","CXCL8 and migration\nshow partial decline",BLUE),("Transient response","Interferon decreases\nmost clearly",TEAL),("Persistent remodeling","LIF/OSM and myeloid\nproteins remain elevated",ORANGE)]
    for i,(h,t,c) in enumerate(boxes):
        y0=.72-i*.31; ax.add_patch(FancyBboxPatch((.05,y0),.9,.23,transform=ax.transAxes,boxstyle="round,pad=.025,rounding_size=.04",fc=c+"14",ec=c,lw=1.2)); ax.text(.1,y0+.15,h,transform=ax.transAxes,fontweight="bold",color=c,fontsize=9); ax.text(.1,y0+.055,t,transform=ax.transAxes,fontsize=8,color=INK)
    fig.suptitle("Longitudinal proteomics resolves distinct recovery trajectories",x=.02,y=.995,ha="left",fontsize=15,fontweight="bold",color=NAVY)
    save(fig,"Figure_3_longitudinal_proteomics")

def figure4():
    expr=read("CIDP Expression Atlas"); effects=read("CIDP Module Effects"); fr=read("CIDP Cell Fractions")
    genes=["CLDN5","ICAM1","VCAM1","LIFR","IL6ST","OSMR","CXCL8","C3","C3AR1","FCGR2A","JUN","SOX10","MPZ"]
    groups=["Macrophage","Granulocyte","BNB_EC","Pericyte","Perineurium","Myelinating_SC","Nonmyelinating_SC","Repair_damage_SC"]
    de=expr[expr.disease.eq("CIDP")&expr.gene.isin(genes)&expr.cell_group.isin(groups)].copy(); write(de,"Figure_4A_CIDP_expression_dotplot")
    mat=de.pivot(index="gene",columns="cell_group",values="mean_log2cpm").reindex(index=genes,columns=groups)
    fig=plt.figure(figsize=(12.6,9.8)); gs=fig.add_gridspec(2,2,width_ratios=[1.0,1.14],height_ratios=[1,1.38],wspace=.78,hspace=.48)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"A","Cell-state localization in CIDP nerve")
    for i,g in enumerate(genes):
        for j,cg in enumerate(groups):
            v=mat.loc[g,cg] if pd.notna(mat.loc[g,cg]) else 0; ax.scatter(j,i,s=max(8,v*8),c=v,cmap="viridis",vmin=0,vmax=max(10,np.nanmax(mat.values)),edgecolor="white",lw=.3)
    ax.set_xticks(range(len(groups)),["Macro","Gran","BNB EC","Pericyte","Perineurium","mySC","nmSC","repairSC"],rotation=35,ha="right"); ax.set_yticks(range(len(genes)),genes); ax.invert_yaxis(); ax.set_xlim(-.6,len(groups)-.4)
    e=effects[effects.comparison.eq("CIDP_vs_CIAP")].copy(); write(e,"Figure_4B_CIDP_module_effects")
    panels=["BNB_identity_integrity","Leukocyte_transmigration","CXCL8_CXCR1_2","LIF_LIFR","Complement","Fc_receptor","Macrophage_state","Schwann_myelin_repair"]
    cgs=["Macrophage","BNB_EC","Perineurium","Myelinating_SC","Nonmyelinating_SC","Repair_damage_SC"]
    pv=e.pivot(index="panel",columns="cell_group",values="delta_mean_z").reindex(panels)[cgs]
    fq=e.pivot(index="panel",columns="cell_group",values="fdr_within_celltype_comparison").reindex(panels)[cgs]
    ax=fig.add_subplot(gs[0,1]); panel(ax,"B","CIDP-versus-CIAP program shifts")
    im=ax.imshow(pv,cmap=mpl.colors.LinearSegmentedColormap.from_list("div",["#2C7BB6","#F7F7F7","#D7191C"]),vmin=-1.2,vmax=1.2,aspect="auto")
    ax.set_xticks(range(len(cgs)),["Macro","BNB EC","Perineurium","mySC","nmSC","repairSC"],rotation=30,ha="right"); ax.set_yticks(range(len(panels)),["BNB identity","Migration","CXCL8","LIF/LIFR","Complement","Fc receptor","Macrophage state","Schwann repair"])
    for i in range(len(panels)):
        for j in range(len(cgs)):
            if pd.notna(fq.iloc[i,j]) and fq.iloc[i,j]<.05: ax.scatter(j+.32,i-.32,marker="D",s=16,color=INK,edgecolor="white",lw=.3)
    cb=fig.colorbar(im,ax=ax,fraction=.036,pad=.025); cb.set_label("Δ module score",fontsize=8); ax.text(.99,-.24,"◆ FDR < 0.05",transform=ax.transAxes,ha="right",fontsize=7,color=GREY)
    f=fr[fr.comparison.eq("CIDP_vs_CIAP")].copy(); write(f,"Figure_4C_CIDP_cell_fractions")
    ax=fig.add_subplot(gs[1,0]); panel(ax,"C","Cellular composition")
    f=f.sort_values("difference_fraction"); y=np.arange(len(f)); ax.axvline(0,color="#9CA3AF",lw=.8); ax.scatter(f.difference_fraction,y,s=34,c=[ORANGE if q<.05 else BLUE for q in f.fdr_across_cell_groups],edgecolor="white",lw=.5)
    ax.set_yticks(y,[x.replace("_"," ") for x in f.cell_group]); ax.set_xlabel("Difference in mean fraction (CIDP − CIAP)"); clean(ax)
    sm=read("CIDP Sample Modules"); cidp=sm[sm.disease.eq("CIDP")&sm.incat.notna()].copy(); sel=cidp[cidp.cell_group.isin(["Macrophage","BNB_EC","Repair_damage_SC"])&cidp.panel.isin(["Complement","Fc_receptor","Leukocyte_transmigration","LIF_LIFR","Schwann_myelin_repair"])]; cor=[]
    for (cg,pn),g in sel.groupby(["cell_group","panel"]):
        if len(g)>=5:
            rho,p=stats.spearmanr(g.score_z,g.incat); cor.append([cg,pn,len(g),rho,p])
    cor=pd.DataFrame(cor,columns=["cell_group","panel","n","rho_INCAT","p_value"]); write(cor,"Figure_4D_CIDP_clinical_correlations")
    ax=fig.add_subplot(gs[1,1]); panel(ax,"D","Program associations with INCAT")
    if len(cor):
        group_short={"BNB_EC":"EC","Macrophage":"Macro","Repair_damage_SC":"repairSC"}
        panel_short={"Complement":"Complement","Fc_receptor":"Fc receptor","LIF_LIFR":"LIF/LIFR","Leukocyte_transmigration":"Transmigration","Schwann_myelin_repair":"Schwann repair"}
        labels=[f"{group_short[r.cell_group]} · {panel_short[r.panel]}" for _,r in cor.iterrows()]
        yy=np.arange(len(cor))[::-1]; ax.axvline(0,color="#9CA3AF",lw=.8); ax.scatter(cor.rho_INCAT,yy,s=40,c=[ORANGE if p<.05 else TEAL for p in cor.p_value],edgecolor="white",lw=.5); ax.set_yticks(yy,labels,fontsize=6.8); ax.set_xlim(-1.05,1.05); ax.set_xlabel("Spearman ρ with INCAT"); clean(ax)
        ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=ORANGE,markeredgecolor="white",label="Nominal P < 0.05"),
                           Line2D([0],[0],marker="o",color="none",markerfacecolor=TEAL,markeredgecolor="white",label="P ≥ 0.05")],
                  frameon=False,fontsize=7.1,loc="upper right",borderaxespad=.2)
    fig.suptitle("CIDP nerve resolves macrophage, BNB endothelial and Schwann-cell programs",x=.02,y=.995,ha="left",fontsize=15,fontweight="bold",color=NAVY)
    save(fig,"Figure_4_CIDP_target_organ")

def figure5():
    cc=read("Cross Compartment"); write(cc,"Figure_5A_cross_compartment_axes")
    sm=read("CIDP Sample Modules"); cidp=sm[sm.disease.eq("CIDP")].copy()
    chosen=[("Macrophage","Complement"),("Macrophage","Fc_receptor"),("BNB_EC","Leukocyte_transmigration"),("BNB_EC","LIF_LIFR"),("Repair_damage_SC","Schwann_myelin_repair")]
    sub=pd.concat([cidp[(cidp.cell_group.eq(cg))&(cidp.panel.eq(pn))] for cg,pn in chosen]); sub["feature"]=[f"{a} · {b}" for a,b in zip(sub.cell_group,sub.panel)]; pv=sub.pivot_table(index="sample",columns="feature",values="score_z"); pv=(pv-pv.mean())/pv.std(ddof=0); write(pv.reset_index(),"Figure_5B_CIDP_donor_heterogeneity")
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,2,width_ratios=[1.12,1],wspace=.42,hspace=.5)
    ax=fig.add_subplot(gs[0,0])
    ax.text(-.11,1.01,"A",transform=ax.transAxes,fontsize=13,fontweight="bold",color=NAVY,va="bottom")
    ax.text(-.015,1.01,"Cross-compartment signaling axes",transform=ax.transAxes,fontsize=10.5,fontweight="bold",color=NAVY,va="bottom")
    axes=cc.axis.tolist(); x=np.arange(len(axes)); ax.bar(x-.18,cc.acute_GBS_hedges_g,width=.36,color=BLUE,label="Acute GBS ligand effect")
    scaled=(cc.CIDP_mean_receptor_log2CPM-cc.CIDP_mean_receptor_log2CPM.mean())/cc.CIDP_mean_receptor_log2CPM.std(ddof=0); ax.bar(x+.18,scaled,width=.36,color=ORANGE,label="CIDP receptor expression (scaled)")
    ax.axhline(0,color="#9CA3AF",lw=.8); ax.set_xticks(x,axes,rotation=28,ha="right"); ax.set_ylabel("Standardized signal"); ax.legend(frameon=False,fontsize=7.5); clean(ax)
    ax=fig.add_subplot(gs[0,1])
    ax.text(-.11,1.01,"B",transform=ax.transAxes,fontsize=13,fontweight="bold",color=NAVY,va="bottom")
    ax.text(-.015,1.01,"CIDP donor heterogeneity",transform=ax.transAxes,fontsize=10.5,fontweight="bold",color=NAVY,va="bottom")
    im=ax.imshow(pv.values,cmap=mpl.colors.LinearSegmentedColormap.from_list("div",["#2C7BB6","#F7F7F7","#D7191C"]),vmin=-2,vmax=2,aspect="auto")
    ax.set_xticks(range(pv.shape[1]),[x.replace("_"," ") for x in pv.columns],rotation=35,ha="right"); ax.set_yticks(range(pv.shape[0]),pv.index,fontsize=6.5); cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.025); cb.set_label("Within-feature z",fontsize=8)
    ax=fig.add_subplot(gs[1,:]);
    ax.text(-.055,.965,"C",transform=ax.transAxes,fontsize=13,fontweight="bold",color=NAVY,va="bottom")
    ax.text(.015,.965,"Myeloid–BNB–Schwann communication model",transform=ax.transAxes,fontsize=10.5,fontweight="bold",color=NAVY,va="bottom")
    ax.set_xlim(0,12); ax.set_ylim(0,3.2); ax.axis("off")
    node_y=.76; node_h=1.05; node_center=node_y+node_h/2
    nodes=[("Circulating\nmyeloid cells",1.2,BLUE),("BNB\nendothelium",4.3,TEAL),("Nerve\nmacrophages",7.4,ORANGE),("Schwann-cell\ninjury and repair",10.6,PURPLE)]
    for t,x0,c in nodes:
        ax.add_patch(FancyBboxPatch((x0-.85,node_y),1.7,node_h,boxstyle="round,pad=.05,rounding_size=.12",fc=c+"18",ec=c,lw=1.5)); ax.text(x0,node_center,t,ha="center",va="center",fontweight="bold",color=c,fontsize=9)
    edges=[(2.15,3.32,"CXCL8 · CCL2\nICAM1/VCAM1"),(5.25,6.42,"transmigration\ncomplement · Fc"),(8.35,9.62,"TNF · IL1B\nCSF1 · LIF/OSM")]
    for x1,x2,t in edges:
        ax.add_patch(FancyArrowPatch((x1,node_center),(x2,node_center),arrowstyle="-|>",mutation_scale=14,lw=1.4,color=GREY)); ax.text((x1+x2)/2,1.54,t,ha="center",va="bottom",fontsize=7.5,color=INK,linespacing=1.12)
    fig.suptitle("Cell-state coupling organizes the CIDP target-organ circuit",x=.02,y=.995,ha="left",fontsize=15,fontweight="bold",color=NAVY)
    save(fig,"Figure_5_cell_communication")

def figure6():
    phase=read("Phase Evidence"); write(phase,"Figure_6A_phase_evidence")
    fig,ax=plt.subplots(figsize=(12.5,6.5)); ax.set_xlim(0,16.5); ax.set_ylim(0,7.2); ax.axis("off")
    ax.text(.35,6.82,"GBS–CIDP immune continuum: recruitment, recovery and target-organ remodeling",fontsize=15,fontweight="bold",color=NAVY)
    centers=[3.45,7.35,11.25,15.15]; colors=[BLUE,TEAL,ORANGE,PURPLE]
    titles=["Acute GBS","Post-acute /\nrecovery","Active chronic\nnerve","Stable / treated\nCIDP"]
    header_y=5.35; header_w=2.35; header_h=.95
    for title,x,c in zip(titles,centers,colors):
        ax.add_patch(FancyBboxPatch((x-header_w/2,header_y-header_h/2),header_w,header_h,boxstyle="round,pad=.06,rounding_size=.1",fc=c+"18",ec=c,lw=1.5))
        ax.text(x,header_y,title,ha="center",va="center",fontweight="bold",color=c,fontsize=9.3,linespacing=1.15)
    arrow_y=header_y
    for x1,x2 in zip(centers[:-1],centers[1:]):
        ax.add_patch(FancyArrowPatch((x1+header_w/2+.12,arrow_y),(x2-header_w/2-.12,arrow_y),arrowstyle="-|>",mutation_scale=14,lw=1.35,color=GREY))
    rows=[("Dominant\ncompartment",["Circulating\nmyeloid cells","Blood and\nserum","Peripheral nerve\nendoneurium","Peripheral nerve\nimmune interface"]),
          ("Leading\nprograms",["CXCL8 recruitment\nComplement and Fc","Interferon decline\nPersistent LIF/OSM","Macrophage remodeling\nBNB antigen response","Treatment-shaped Fc\nRepair programs"]),
          ("Clinical\ncontext",["Rapid weakness","Recovery trajectory","Activity or relapse","IVIg-treated state"])]
    yvals=[3.82,2.48,1.15]; cell_w=2.62; cell_h=.92
    for (row,vals),y in zip(rows,yvals):
        ax.text(1.86,y,row,fontweight="bold",color=NAVY,ha="right",va="center",fontsize=8.3,linespacing=1.15)
        for x,c,val in zip(centers,colors,vals):
            ax.add_patch(FancyBboxPatch((x-cell_w/2,y-cell_h/2),cell_w,cell_h,boxstyle="round,pad=.04,rounding_size=.06",fc="white",ec=c+"55",lw=1.05))
            ax.text(x,y,val,ha="center",va="center",fontsize=7.7,color=INK,linespacing=1.18)
    save(fig,"Figure_6_integrated_mechanism")

def tables():
    read("Dataset Design").to_csv(RAW/"Table_1_datasets_and_design.csv",index=False)
    read("Cross Compartment").to_csv(RAW/"Table_2_cross_compartment_axes.csv",index=False)
    meta=pd.read_csv(RAW/"Figure_2B_random_effects_meta.csv")
    key=meta[["module","effect","low","high","p","I2","k"]].copy(); key.columns=["program","summary_hedges_g","ci_low","ci_high","hartung_knapp_p","I2_percent","n_cohorts"]
    key.to_csv(RAW/"Table_3_crosscohort_robustness.csv",index=False)
    manifest=[]
    for p in sorted(RAW.glob("*.csv")):
        d=pd.read_csv(p); manifest.append(dict(file=p.name,rows=len(d),columns=len(d.columns),purpose=p.stem.replace("_"," ")))
    pd.DataFrame(manifest).to_csv(RAW/"source_data_manifest.csv",index=False)

def main():
    import sys
    jobs={"1":figure1,"2":figure2,"3":figure3,"4":figure4,"5":figure5,"6":figure6}
    selected=sys.argv[1:] or list(jobs)
    for key in selected:
        print(f"Rendering Figure {key}...",flush=True)
        jobs[key]()
    tables()
    print(json.dumps({"figures":len(list(OUT.glob('*.png'))),"source_tables":len(list(RAW.glob('*.csv')))},indent=2))

if __name__ == "__main__": main()
