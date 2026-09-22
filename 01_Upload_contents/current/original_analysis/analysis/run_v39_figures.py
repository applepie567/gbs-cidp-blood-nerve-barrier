"""Reproduce the revision-specific CSF correction and all seven figures.

Run from any directory with the original project dependencies available.
This does not reprocess original matrices or recompute donor-level nerve analyses.
"""
from pathlib import Path
import importlib.util
import subprocess
import sys

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('v39_validation',HERE/'08_extended_validation.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
result=module.csf_meta()
assert len(result)==10
assert set(result.protein)=={'SELE','IL2RA','CCL3','CR2','CD1C','THBD','NRP1','CD38','IL6','CD5'}
for name in ['10_figures_1_2_5.py','11_figure_3.py','12_figure_4.py','14_supplementary_figures_v39.py']:
    subprocess.run([sys.executable,str(HERE/name)],check=True)
print('v39 CSF correction and seven figures complete')
