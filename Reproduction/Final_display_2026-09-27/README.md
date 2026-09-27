# Final figure display revision

Requirements: Python, matplotlib, numpy and pandas.

From this directory, run:

```bash
python render_final_figures.py --output ../../Figures
```

This rebuilds Figure 1, Figure 5 and Supplementary Figures S5–S7 from the supplied study catalogue and frozen aggregate CSVs. It changes layout, keys and redundant text only. Geometry assertions verify that display edits do not change numerical plot coordinates or heatmap values. The base script retains historical figure numbering; the wrapper maps its Figure 4 to current Figure 5. All other figures retain their existing exports and source scripts, as recorded in `Source_data/Figure_source_index.tsv`.
