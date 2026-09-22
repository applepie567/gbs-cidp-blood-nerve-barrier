"""Re-export the per-panel figure tables from included upstream analysis outputs.

Run from any working directory. This preserves the fixed manuscript design;
change the design and table definitions when a future analysis changes them.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
C=Path(__file__).resolve().parents[1]
def rd(p): return pd.read_csv(p)
def rel(p): return str(p.relative_to(C))
def js(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
L=C/'legacy_reference/analysis_update/source_data';F=C/'legacy_reference/strengthening_v41/results';S=C/'results'
N=C/'followup_20260921/results/nerve_multiplicity_91.csv'
idx=[];sheets=[]
def add(name,figure,panel,df,source,unit,note,kind='Aggregate analysis estimate',url=''):
    df=df.reset_index(drop=True)
    first=['display_model','analysis','model'] if name=='F4D_Fc_models' else ['cell_group','program','model'] if name in ['F5A_Endothelial','F5B_Schwann'] else []
    df=df[[x for x in first if x in df]+[x for x in df if x not in first]]
    df.to_csv(C/'figure_source_data'/f'{name}.csv',index=False)
    idx.append(dict(sheet=name,figure=figure,panel=panel,rows=len(df),csv=f'{name}.csv',data_level=kind,unit=unit,upstream_file=source,source_url=url,notes=note))
    payload=json.loads(df.to_json(orient='split',index=False));sheets.append(dict(name=name,columns=payload['columns'],rows=payload['data']))
nerve='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285983'
blood='GSE211225; GSE31014; PRJNA1293757'
csf='https://doi.org/10.3389/fimmu.2023.1241199'
add('F1_Study_design','Figure 1','All',pd.DataFrame([
 {'compartment':'Blood','comparison':'Acute GBS versus healthy controls','cohorts':3,'participants':31,'programme_tests':7,'source':blood},
 {'compartment':'CSF','comparison':'GBS versus HC; GBS versus CIDP; CIDP versus HC','cohorts':2,'participants':None,'paired_proteins':11,'additional_comparison_estimates':16,'source':csf},
 {'compartment':'Peripheral nerve','comparison':'CIDP versus CIAP','cohorts':1,'participants':20,'CIDP':9,'CIAP':11,'source':nerve}]),'Manuscript; analysis/make_figures.py','Study design and counts','Schematic, not experimental measurements. CSF contrasts share participants.','Study design')
for nm,panel,fn,unit,note in [
 ('F2A_Blood_cohort','A','Blood_cohort_effects.csv','Hedges g','Seven programmes in three cohorts; FDR retained from each cohort.'),
 ('F2B_Blood_meta','B','Blood_meta.csv','Hedges g','DerSimonian–Laird pooled effects with modified Hartung–Knapp intervals.'),
 ('F2C_Blood_leaveout','C','Figure_2C_leave_one_cohort_out.csv','Hedges g','Pooled effects after each omitted cohort; these are sensitivity estimates, not additional cohorts.')]:
    add(nm,'Figure 2',panel,rd(L/fn),rel(L/fn),unit,note,url=blood)
m=rd(S/'CSF_meta_11_proteins.csv');order=['IL8','SELE','IL2RA','CCL3','CR2','CD1C','THBD','NRP1','CD38','IL6','CD5'];m=m.set_index('protein').loc[order].reset_index()
rows=[];pooled=[]
for _,r in m.iterrows():
    for cohort in ['discovery','replication']:
        rows.append(dict(protein=r.protein,cohort=cohort,beta=r[cohort+'_beta'],se=r[cohort+'_se'],ci_low=r[cohort+'_beta']-1.95996398454*r[cohort+'_se'],ci_high=r[cohort+'_beta']+1.95996398454*r[cohort+'_se'],source_p=r[cohort+'_source_p'],source_q=r[cohort+'_source_q'],source_table=r.source_table))
    for model in ['fixed','random']:
        pooled.append(dict(protein=r.protein,model=model,beta=r[model+'_beta'],se=r[model+'_se'],ci_low=r[model+'_low'],ci_high=r[model+'_high'],p_value=r[model+'_p'],q_BH_11=r['q_'+model+'_BH_11']))
add('F3A_CSF_cohorts','Figure 3','A',pd.DataFrame(rows),rel(S/'CSF_meta_11_proteins.csv'),'Adjusted log2 NPX coefficient','Normal-approximation display intervals. Source P and q are kept separately.','Published regression estimate',csf)
add('F3B_CSF_pooled','Figure 3','B',pd.DataFrame(pooled),rel(S/'CSF_meta_11_proteins.csv'),'Adjusted log2 NPX coefficient','Plotted fixed and random effects; normal-approximation intervals.',url=csf)
add('F3_Pooling_sensitivity','Figure 3','B and legend',m,rel(S/'CSF_meta_11_proteins.csv'),'Adjusted log2 NPX coefficient','All 11 modified Hartung–Knapp intervals include zero; retain all pooling methods.',url=csf)
d=rd(S/'CSF_disease_comparisons.csv')
for nm,panel,contrast in [('F3C_GBS_CIDP','C','GBS versus CIDP'),('F3D_CIDP_HC','D','CIDP versus HC')]:
    t=d[d.contrast.eq(contrast)].copy();t['orange_source_q_below_05']=t.source_q.lt(.05)
    add(nm,'Figure 3',panel,t,rel(S/'CSF_disease_comparisons.csv'),'Adjusted log2 NPX coefficient','Intervals from normal approximation; colours reflect source q. Shared participants.','Published regression estimate',csf)
add('F4A_Nerve_expression','Figure 4','A',rd(L/'Figure_4A_mean_expression.csv'),rel(L/'Figure_4A_mean_expression.csv'),'Mean donor log2(CPM + 0.5)','CIDP; thirteen genes in seven cell groups.',url=nerve)
d=rd(L/'CIDP_module_effects.csv');q91=rd(N)
d=d.merge(q91[['comparison','cell_group','panel','q_bh_91']],on=['comparison','cell_group','panel'],validate='one_to_one')
cells=['Macrophage','BNB_EC','Perineurium','Myelinating_SC','Nonmyelinating_SC','Repair_damage_SC']
d=d[d.comparison.eq('CIDP_vs_CIAP')&d.cell_group.isin(cells)]
add('F4B_Module_contrasts','Figure 4','B',d,rel(L/'CIDP_module_effects.csv')+'; '+rel(N),'Mean module z-score difference','Diamonds use BH within cell group (eight modules). Expanded q91 supplied separately.',url=nerve)
add('F4C_Cell_fractions','Figure 4','C',rd(L/'CIDP_cell_fractions.csv'),rel(L/'CIDP_cell_fractions.csv'),'Fraction difference, CIDP minus CIAP','Relative representation among sampled nuclei; not absolute tissue cell abundance.',url=nerve)
sens=rd(F/'fc_sensitivity.csv');comp=rd(F/'fc_composition_adjustment.csv')
d=pd.concat([sens[sens.analysis.isin(['All donors','Age and center'])],comp],ignore_index=True)
d['display_model']=['Unadjusted','Age and centre','Age, centre and composition']
add('F4D_Fc_models','Figure 4','D',d,rel(F/'fc_sensitivity.csv')+'; '+rel(F/'fc_composition_adjustment.csv'),'Mean module z-score difference','HC3 intervals. Adjustment for state composition attenuates the association.',url=nerve)
u=rd(S/'nerve_program_contrasts.csv');a=rd(S/'nerve_program_adjusted.csv')
u=u.merge(q91[q91.source_family.eq('Nine functional programme extensions')][['cell_group','panel','q_bh_91']].rename(columns={'panel':'program'}),on=['cell_group','program'],validate='one_to_one');u['model']='Unadjusted'
d=pd.concat([u,a],ignore_index=True)
for name,panel,mask in [('F5A_Endothelial','A',d.cell_group.eq('BNB_EC')),('F5B_Schwann','B',~d.cell_group.eq('BNB_EC'))]:
    add(name,'Figure 5',panel,d[mask],rel(S/'nerve_program_contrasts.csv')+'; '+rel(S/'nerve_program_adjusted.csv')+'; '+rel(N),'Mean programme z-score difference','Unadjusted primary P: Mann–Whitney. Display intervals: HC3. Adjusted P/q: HC3 family of nine. q89 is historical; q91 is the current expanded family.',url=nerve)
coords=rd(C/'source_data/nerve_clinical_plot_data.csv')
for name,panel,prog in [('F5C_Myelin_points','C','Myelin maintenance'),('F5D_Injury_points','D','Injury response')]:
    add(name,'Figure 5',panel,coords[coords.program.eq(prog)],'source_data/nerve_clinical_plot_data.csv','Programme z-score; INCAT disability score','Nine CIDP donors in pooled repair/damage Schwann cells. Row order does not identify or link donors across programmes.','Anonymous donor-level plot value',nerve)
add('F5CD_Clinical_statistics','Figure 5','C and D',rd(S/'nerve_clinical_correlations.csv'),rel(S/'nerve_clinical_correlations.csv'),'Spearman rho; P and q','Exact P enumerates all 9! assignments. Expanded q17 combines 15 original asymptotic P values with two new exact P values.',url=nerve)
f=C/'original_analysis/results/tables/genetic_cell_localization_bootstrap.csv'
add('S1_Bootstrap_expression','Supplementary Figure 1','All',rd(f),rel(f),'Mean donor log2(CP10k + 1); percent expressing','Donor-bootstrap expression intervals in 15 cell groups for four genes. Per-nucleus workflow without module eligibility threshold.',url=nerve)
url='https://doi.org/10.1007/s00401-025-02936-w';f=C/'original_analysis/results/tables/cidp_external_nerve_validation.csv';d=rd(f)
t=d.iloc[[0]].copy();t['positive_biopsies']=52;t['total_biopsies']=55
add('S2_Complement_count','Supplementary Figure 2','Proportion',t,rel(f),'Proportion (0 to 1); displayed as percent','Published count; Wilson 95% interval. No individual biopsy observations are reconstructed.','Published count',url)
add('S2_Tissue_observations','Supplementary Figure 2','Text',d,rel(f),'Published observations','Different endpoints from macrophage Fc expression; biological context, not replication of the Fc association.','Published observations',url)
add('S3A_Fc_genes','Supplementary Figure 3','A',rd(F/'fc_component_genes.csv'),rel(F/'fc_component_genes.csv'),'Difference in mean log2(CPM + 0.5)','HC3 intervals; q across eleven available genes.',url=nerve)
co=comp.copy();co['analysis']='Age, center and composition';d=pd.concat([sens.iloc[:5],co],ignore_index=True)
add('S3B_Fc_sensitivity','Supplementary Figure 3','B',d,rel(F/'fc_sensitivity.csv')+'; '+rel(F/'fc_composition_adjustment.csv'),'Mean module z-score difference','Six prespecified displayed rows; donor/model sensitivities within the same cohort.',url=nerve)
add('S4_Fc_annotations','Supplementary Figure 4','All',rd(F/'fc_annotation_sensitivity.csv'),rel(F/'fc_annotation_sensitivity.csv'),'Mean module z-score difference','Three broad definitions and nine eligible refined subtypes; donor counts vary by annotation.',url=nerve)
add('S5A_Gene_localisation','Supplementary Figure 5','A',rd(L/'Figure_5B_donor_mean_expression_and_detection.csv'),rel(L/'Figure_5B_donor_mean_expression_and_detection.csv'),'Mean donor log2(CP10k + 1); percent expressing','Dot colour = expression; area = 5 × percent. Zero detection is shown by a cross.',url=nerve)
t=rd(L/'CIDP_genetic_evidence.csv');t['odds_ratio']=np.nan;t['odds_ratio_ci_low']=np.nan;t['odds_ratio_ci_high']=np.nan;t.loc[t.gene.eq('CDH4'),['odds_ratio','odds_ratio_ci_low','odds_ratio_ci_high']]=[4.37,2.61,7.33]
add('S5B_Genetic_evidence','Supplementary Figure 5','B',t,rel(L/'CIDP_genetic_evidence.csv'),'Odds ratio; MR beta; posterior probability','CDH4 association in women. MR estimates are copied from authoritative input rows, including GNG7 −2.11 and SLC39A3 −1.42.','Published association estimate','https://doi.org/10.1016/j.xhgg.2024.100317')
add('S5C_Gene_contrasts','Supplementary Figure 5','C',rd(L/'Figure_5D_genetic_expression_contrasts.csv'),rel(L/'Figure_5D_genetic_expression_contrasts.csv'),'Mean log2(CP10k + 1) difference','CIDP minus CIAP in ten cell groups; no displayed contrast passes correction.',url=nerve)
d=rd(L/'CIDP_clinical.csv');p=d.p_value.to_numpy();o=np.argsort(p);qq=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[o]=np.minimum(qq,1);d['q_BH_15']=q
add('S6_Clinical_correlations','Supplementary Figure 6','All',d,rel(L/'CIDP_clinical.csv'),'Spearman rho','Original 15 asymptotic correlations and BH15; all q > 0.05.',url=nerve)
d=rd(S/'BNB_normal_reference_expression.csv');d['plotted_log2_fpkm_plus_1']=np.log2(d.fpkm+1)
add('S7_Normal_BNB_expression','Supplementary Figure 7','All',d,rel(S/'BNB_normal_reference_expression.csv'),'FPKM; plotted log2(FPKM + 1)','Eleven genes × six preparations. Normal anatomical reference; no disease comparison.','Published expression value','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107574')
add('Support_q91','Supporting data','Multiplicity',q91,rel(N),'P; Benjamini–Hochberg q','80 original module tests + 9 programme extensions + 2 graft-reference contrasts.',url=nerve)
add('Support_gene_coverage','Supporting data','Programmes',rd(S/'nerve_gene_coverage.csv'),rel(S/'nerve_gene_coverage.csv'),'Gene availability and standardisation','Measured, variable genes used in each programme; missing genes are not imputed.',url=nerve)
add('Public_input_manifest','Supporting data','Primary inputs',rd(C/'source_data/public_input_manifest.csv'),'source_data/public_input_manifest.csv','Files, URLs and SHA256','Public sequencing/count and normal BNB sources. Download separately for reanalysis.','Primary input manifest')
index=pd.DataFrame(idx);index.to_csv(C/'figure_source_data/INDEX.csv',index=False)
notes=pd.DataFrame([
 ('Manuscript','Immune recruitment signatures and cellular context in inflammatory neuropathies'),
 ('Prepared version','v2.4.0; 2026-09-22; prepared update, not an already published release'),
 ('Scope','Source data for all five main figures and seven supplementary figures in the current manuscript.'),
 ('What these data are','The exact numerical inputs or reported source estimates used for plotted points, intervals, heatmaps and counts.'),
 ('Primary raw data','FASTQ, CellBender matrices and original proteomics acquisitions remain at their public repositories. They are not recreated in this workbook.'),
 ('Figure 1','A study-design schematic. Its table contains design labels and counts, not experimental measurements.'),
 ('Missing values','Blank cells mean unavailable or not applicable, not zero. Read each sheet with INDEX.'),
 ('Clinical coordinates','Figure 5 C/D contain anonymous score/INCAT pairs; no donor linkage across panels is implied by row order.'),
 ('Current correction families','Preserve distinct q8, q80/q89 historical, q91 expanded, q9 new and q17 clinical families; they are not interchangeable.'),
 ('Protein validation','Public quantitative exports and feasibility audit are included in the GitHub ZIP. Missing reliable sample-to-diagnosis mapping prevents independent nerve validation; this is not a negative validation result.'),
 ('Reproduction','Run python current/analysis/reproduce_figures.py --out rebuilt_figures from the extracted repository root.'),
 ('Data licence','Project-derived tables and figures: CC BY 4.0. External source material retains original terms.')],columns=['topic','details'])
notes.to_csv(C/'figure_source_data/README.csv',index=False)
def split(df,name):
    o=json.loads(df.to_json(orient='split',index=False));return {'name':name,'columns':o['columns'],'rows':o['data']}
# Use Python-native float serialization, not pandas' default 10-digit JSON rounding.
def exact_payload(df,name):
    rows=[[None if pd.isna(v) else v.item() if isinstance(v,np.generic) else v for v in row] for row in df.itertuples(index=False,name=None)]
    return {'name':name,'columns':list(df.columns),'rows':rows}
all_sheets=[exact_payload(notes,'README'),exact_payload(index,'INDEX')]
for x in sheets:all_sheets.append(exact_payload(rd(C/'figure_source_data'/f"{x['name']}.csv"),x['name']))
js(C/'figure_source_data/workbook_payload.json',all_sheets)

print('Exported per-panel CSVs and workbook payload; the frozen workbook is not overwritten.')
