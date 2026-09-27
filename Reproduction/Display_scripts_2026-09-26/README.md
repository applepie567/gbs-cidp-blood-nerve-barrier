# Current figure displays

Run from this directory with the numerical plotting dependencies from the archived environment:

```bash
python render_display.py --output figures
```

This renders Figure 1, Figure 3 and Supplementary Figure S8. It uses `base_make_figures.py`, copied unchanged from the v2.4.0 release, and the six unchanged aggregate CSV files in `data/`. The script checks that legend and colour edits leave plotted coordinates, effect estimates, confidence intervals and axis limits unchanged. No statistical analysis is rerun.

Figure 1 now presents the study questions. Figure 3 distinguishes source q thresholds in panels C and D. Supplementary Figure S8 is the former main Figure 5, with the model legend above A and B and grey points in C and D. Other figures are supplied as exact copies of the manuscript images, with reconstruction information in Supplementary Data 2.

The 26 September 2026 Figure 1 adds GSE213455 (4 CIDP and 9 vasculitic neuropathy biopsies). Figure S9 and its complete analysis are reproduced from Supplementary Data 4 with `python render_external_nerve.py --root . --output figures`. The old plotting inputs and numerical results are unchanged.
