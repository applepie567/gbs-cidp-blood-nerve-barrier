# GitHub 更新说明

本包准备版本号为 v2.4.0，日期为 2026-09-22；尚未上传 GitHub 或发布 Zenodo。

1. 解压本 ZIP，打开其中的 `GitHub_update_v2.4.0` 文件夹。
2. 将该文件夹内的文件和 `current` 文件夹复制到现有仓库根目录，保留层级。顶层 README、CITATION、元数据及依赖文件会更新；本次稿件、数据和图片均位于 `current/`。
3. 旧的 `analysis_update/`、`strengthening_v41/` 和旧稿件可保留作历史记录。新版 README 已指向 `current/`，请从这里取本次投稿文件。
4. 本地运行 `python current/analysis/verify_release.py`，必要时运行 `python current/analysis/reproduce_figures.py --out rebuilt_figures`。
5. 检查提交差异后提交并推送。上传代码时应上传解压内容；ZIP 可另作为 Release 附件。使用 Git 提交文件夹比网页逐个上传方便。
6. 如随后发布 Release，可使用尚未存在的 `v2.4.0` 标签。取得真实的新版本 DOI 后，再更新 Word 稿件 Data availability 和仓库引文；本包没有预填未发布 DOI。

`current/figure_source_data/` 包含 Excel、逐图逐分面 CSV 和数据索引；`current/figures/` 为当前文稿插图。图1、3、5沿用最后一次删除图内说明后的版本。主稿及补充文件 1–3 均保留当前内容。

蛋白组定量导出位于 `current/protein_validation_20260921/results/`。由于缺少可靠的样本与诊断对应关系，这部分仅完成公开数据导出及可用性评估，不能写成独立验证完成或验证阴性。

“绘图源数据”包括实际绘图使用的坐标、汇总估计和已发表数值；原始 FASTQ、计数矩阵和质谱采集文件需按公开来源另行下载。旧补充分析 ZIP 的内容已从完整保留文件重新组装，本更新包不嵌套该损坏压缩包。
