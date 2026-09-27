"""Build the self-contained composition package and update the full manuscript bundle."""
from pathlib import Path
import shutil,gzip,hashlib,json,csv,platform,sys
from zipfile import ZipFile,ZIP_DEFLATED
import numpy,pandas,scipy,matplotlib

R=Path(__file__).resolve().parents[1];O=R/'outputs';D=O/'Supplementary_Data_6'
P=O/'GBS_CIDP_strengthened_2026-09-27'
OLD=R.parent/'validation_20260927/outputs/GBS_CIDP_external_blood_2026-09-27'
D.mkdir(exist_ok=True)
for folder in ['scripts','results','Figures','input/prior/External_nerve_20260926/inputs','input/prior/External_nerve_20260926/results']:(D/folder).mkdir(exist_ok=True,parents=True)
for n in ['ANALYSIS_PLAN.md','source_manifest.csv']:shutil.copy2(R/n,D/n)
for n in ['download_reference.py','build_reference.py','analyse_composition.py','make_figure.py','verify_models.R']:
    shutil.copy2(R/'scripts'/n,D/'scripts'/n)
prior=R/'input/prior/External_nerve_20260926';new=D/'input/prior/External_nerve_20260926'
for n in ['analyse_external_nerve.py','hdf_reader.py','input_manifest.csv']:shutil.copy2(prior/n,new/n)
for n in ['GSE285983_macrophage_Fc_counts.csv','Fc_donor_scores.csv','GSE213455_sample_mapping.csv']:
    shutil.copy2(prior/'results'/n,new/'results'/n)
shutil.copy2(prior/'inputs/GSE213455_series_matrix.txt.gz',new/'inputs/GSE213455_series_matrix.txt.gz')
with (prior/'inputs/GPL13369_annotation.tsv').open('rb') as src,gzip.open(new/'inputs/GPL13369_annotation.tsv.gz','wb') as dst:shutil.copyfileobj(src,dst)
for p in (R/'results').iterdir():
    if p.suffix in ['.csv','.npz','.json'] or p.name.endswith('.csv.gz') or p.name=='base_R_verification.txt':
        shutil.copy2(p,D/'results'/p.name)
for p in (O/'Figures').iterdir():shutil.copy2(p,D/'Figures'/p.name)
(D/'requirements.txt').write_text('numpy\npandas\nscipy\nmatplotlib\nh5py\n')
(D/'software_versions.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=numpy.__version__,pandas=pandas.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),indent=2))
(D/'README.md').write_text('''# Nerve macrophage marker analysis

This supplement extends the GSE213455 whole nerve analysis. It does not provide independent cell-resolved CIDP versus CIAP replication.

## Result

The reference contains 365,708 author-retained nuclei from 37 GSE285983 donors, aggregated into 13 broad lineages. The 20-marker index excludes all 11 Fc outcome genes. Leave-one-donor-out Spearman correlations are 0.933381 with the observed macrophage nucleus fraction and 0.950450 with the macrophage share of captured counts. Leaving out entire centres gives 0.911332 and 0.927217, respectively. All within-centre checks remain positive and the recorded feasibility criteria are met.

In the independent four-CIDP versus nine-vasculitic-neuropathy cohort, the fixed marker index correlates with the Fc score at rho 0.906593. The original difference is -0.836553 (HC3 95% CI -1.763592 to 0.090486). Adding the index gives 0.208146 (-0.915957 to 1.332250, model P 0.688623). Mean-probe and donor-omission sensitivities are supplied in full. The conditional result concerns marker expression and does not separate cell abundance from cell state. The index is not an absolute cell fraction, and same-platform atlas calibration does not validate cross-platform accuracy.

## Reproduce the statistical analysis

From this directory, run:

```
python scripts/analyse_composition.py
python scripts/make_figure.py
Rscript scripts/verify_models.R .
```

The analysis uses the included all-gene donor-lineage counts, metadata, external series matrix, GPL13369 annotation, and unchanged external Fc scores. It needs no large matrix download. Gene ranks use a fixed common background of 16,633 genes, excluding the 11 outcome genes. Marker selection is rebuilt within every held-out training set. `results/all_marker_selections.csv` retains every selected gene and selection statistic. The NPZ array has rows matching `atlas_donor_lineage_metadata.csv`, and its `genes` entry identifies columns. Counts include all deposited gene symbols after summing duplicates. The original 11-gene Fc reconstruction is exact for every donor.

To repeat aggregation from the public source matrices, install h5py or use an existing compatible HDF5 library and run:

```
python scripts/download_reference.py
python scripts/build_reference.py
python scripts/analyse_composition.py
```

The download script retrieves and verifies the 37 matrices and full original nucleus metadata. Source URLs and hashes are in `source_manifest.csv`. The compact `atlas_nucleus_annotation.csv.gz` records the retained nucleus labels used here. The full original metadata and large HDF5 files are not redistributed in this package.

## Figure and table sources

- S12A: `results/donor_holdout_predictions.csv`
- S12B: `results/centre_holdout_predictions.csv`
- S12C: `results/external_donor_scores.csv`
- S12D and Table S6: `results/external_conditional_models.csv`
- Additional calibration metrics: `results/feasibility_metrics.csv`
- Marker identities: `results/frozen_macrophage_markers.csv`
- Omission sensitivity: `results/external_leave_one_out.csv`

All points are donors, not individual nuclei. The original Fc score is unchanged. Conditional models use the existing cohort score scale, including during omission checks. Model intervals and P values are separate from the original exact rank P value.

## Provenance and availability

Heming et al., Nature Communications 2025, doi:10.1038/s41467-025-62964-8, GSE285983.
Cerri et al., Frontiers in Oncology 2022, doi:10.3389/fonc.2022.974751, GSE213455 and GPL13369.

This is an exploratory extension specified after the original external disease result was known. Its marker and feasibility rules were written before deriving the new marker results. Public data retain their source attribution and terms. The previous public Zenodo v2.4.0 archive does not contain this analysis. No new DOI is claimed.
''')
shutil.copy2(R.parent/'validation_20260927/outputs/Supplementary_Data_5/LICENSE',D/'LICENSE')

def archive(src,dest):
    with ZipFile(dest,'w',ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(src.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:z.write(p,src.name+'/'+str(p.relative_to(src)))
    with ZipFile(dest) as z:assert z.testzip() is None
    print(dest.name,dest.stat().st_size,flush=True)

def manifest(root):
    file=root/'SHA256SUMS.txt'
    lines=[]
    for p in sorted(root.rglob('*')):
        if p.is_file() and p!=file and '__pycache__' not in p.parts:
            with p.open('rb') as f: h=hashlib.file_digest(f,'sha256').hexdigest()
            lines.append(h+'  '+str(p.relative_to(root)))
    file.write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
    manifest(D);archive(D,O/'Supplementary_Data_6.zip')
    shutil.copytree(OLD,P,dirs_exist_ok=True)
    for n in ['Manuscript_refocused.docx','Supplementary_Information_refocused.docx','Supplementary_Data_6.zip']:shutil.copy2(O/n,P/n)
    for p in (O/'Figures').iterdir():shutil.copy2(p,P/'Figures'/p.name)
    source=P/'Source_data/Nerve_composition';source.mkdir(exist_ok=True)
    for n in ['donor_holdout_predictions.csv','centre_holdout_predictions.csv','external_donor_scores.csv','external_conditional_models.csv','feasibility_metrics.csv','frozen_macrophage_markers.csv','external_leave_one_out.csv']:
        shutil.copy2(R/'results'/n,source/n)
    idx=P/'Source_data/Figure_source_index.tsv'
    with idx.open() as f:reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
    for panel,name in [('S12A','donor_holdout_predictions.csv'),('S12B','centre_holdout_predictions.csv'),('S12C','external_donor_scores.csv'),('S12D','external_conditional_models.csv')]:
        rows.append(dict(zip(fields,[f'Supplementary Figure {panel}','Not in v2.4.0','Figures/Supplementary_Figure_12.png','Source_data/Nerve_composition/'+name,'New donor and centre calibration or external composition sensitivity'])))
    rows.append(dict(zip(fields,['Supplementary Table S6','Not in v2.4.0','Supplementary_Information_refocused.docx','Source_data/Nerve_composition/external_conditional_models.csv','Conditional models use the original Fc scale'])))
    with idx.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    old_readme=(P/'README.md').read_text()
    (P/'README_external_blood_2026-09-27.md').write_text(old_readme)
    (P/'README.md').write_text('''# Manuscript update with nerve composition analysis

Use the two Word files at this directory level as the current manuscript and supplement. Changed text is highlighted yellow. There are five main figures, twelve supplementary figures, two main tables and six supplementary tables. The new analysis is Figure S12 and Table S6. All preceding blood and nerve results are retained.

Supplementary Data 6 contains the new reference counts, held-out marker sets, calibration, external models and executable code. Its README gives the reproduction commands. Supplementary Data 1 to 5 retain their previous scope and numerical contents. All current figures are in `Figures/`, and every new panel has source CSVs in `Source_data/Nerve_composition/`. The figure-source index maps current and historical labels.

The macrophage marker index tracks observed atlas composition in held-out donors. In the independent nerve biopsies, adjustment for this index changes the Fc difference from -0.837 to 0.208 (95% CI -0.916 to 1.332). This strengthens interpretation of the tissue result but does not supply independent cell-resolved CIDP versus CIAP replication. The index may reflect both abundance and cell state. The detailed Chinese change note explains the result.

The previous external blood update is documented in `README_external_blood_2026-09-27.md` and `CHANGES_2026-09-27_zh.md`. Its references to eleven supplementary figures and five supplementary tables describe the preceding version. Historical archives retain older text and figures for provenance.

These new analyses have not been published as a new GitHub or Zenodo release. DOI https://doi.org/10.5281/zenodo.22885750 remains the frozen v2.4.0 archive.
''')
    (P/'CHANGES_nerve_composition_2026-09-27_zh.md').write_text('''# 神经巨噬细胞组成补强

本轮完成了实际分析，并已写入正文和补充材料。修改处标黄。

共重新汇总37名供者的365708个细胞核，形成13类细胞的供者水平全基因计数。另建20个巨噬细胞标志基因的表达指数，排除了原Fc模块的全部11个基因。每次留出整名供者，重新选择标志基因。指数与观察到的巨噬细胞核比例相关系数为0.933，与其计数占比为0.950。按中心留出后分别为0.911和0.927，符合事先记录的可行性条件。

在独立GSE213455队列的4例CIDP和9例血管炎性神经病中，指数与Fc评分的相关系数为0.907。原Fc差异为−0.837，加入该指数后的差异为0.208，95%置信区间为−0.916至1.332，模型P值为0.689。平均探针和逐一剔除供者的结果均完整保留。

这说明独立组织队列的Fc差异对更广泛的巨噬细胞表达贡献敏感。该指数不是细胞比例的直接测量，其变化也可能包含细胞状态变化。因此，本轮增加了对组织信号来源的解释，尚未完成与原CIDP对CIAP比较匹配的独立巨噬细胞验证。

新增补充图S12、三线表S6和Supplementary Data 6。主要系数、HC3置信区间及P值已通过独立R实现复核。正文的摘要、结果、讨论、方法及数据代码声明已同步更新。原血液、CSF及其他神经结果保留。新的图、源数据和分析代码均在本包内。

请使用最外层两个Word文件作为当前稿件。此次内容尚未发布为新的GitHub或Zenodo版本，现有v2.4.0 DOI不包含本轮新结果。
''')
    # Retain the scripts that assembled the reviewed documents without mixing them with the scientific pipeline.
    assembly=P/'Reproduction/Composition_document_update';assembly.mkdir(exist_ok=True)
    for n in ['update_manuscripts.py','package_results.py']:shutil.copy2(R/'scripts'/n,assembly/n)
    (assembly/'README.md').write_text('These scripts record the document and package assembly in the working directory used for this revision. The portable scientific reproduction scripts and all numerical inputs are in Supplementary Data 6.\n')
    manifest(P);archive(P,O/'GBS_CIDP_strengthened_2026-09-27.zip')
