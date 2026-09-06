# Release checks

Run from the repository root:

```bash
python analysis/run_v2_release.py
python analysis/render_current_figures.py
python tests/validate_release.py
python tests/validate_public_distribution.py
```

The first command checks aggregate-supported calculations and workbook/figure alignment. Rendering regenerates five main figures. The next validation compares all seven figure PNGs with manuscript v37 and checks all 30 workbook sheets against their CSV sources. The distribution check requires omitted individual-level intermediates to be absent.

Full donor recalculation is available separately with `--with-local-donor-tables` after reconstructing required local inputs. It is not part of public aggregate verification and was not run for this update.

`metadata/SHA256SUMS.txt` lists distributed files. These hashes describe the release package, not the original database downloads. Historical v36 reports are in `metadata/history/`; current reports state their own computational scope.
