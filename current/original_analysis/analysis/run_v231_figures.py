"""Regenerate revised Figures 3 to 5 from public aggregate inputs."""
from pathlib import Path
import subprocess, sys
root=Path(__file__).resolve().parent
for name in ['04_compile_published_csf_evidence.py','06_build_cross_compartment_tables.py','11_figure_3.py','redraw_figures_4_5.py']:
    subprocess.run([sys.executable,str(root/name)],check=True)
