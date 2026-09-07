"""Regenerate the five v36 main figures from the supplied inputs."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name in ['10_figures_1_2_5.py','11_figure_3.py','12_figure_4.py']:
    subprocess.run([sys.executable,str(root/name)],check=True)
