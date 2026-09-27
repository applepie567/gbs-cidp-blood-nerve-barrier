"""Render the current study diagram and colour legends from unchanged aggregate data.

Run beside base_make_figures.py and data/. No statistical analyses are rerun.
"""
from pathlib import Path
import importlib.util
import argparse
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Bbox

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('base_figures',ROOT/'base_make_figures.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)

def geometry(fig):
    out=[]
    for ax in fig.axes:
        out.extend([repr(ax.get_xlim()),repr(ax.get_ylim())])
        for line in ax.lines:
            out.extend([np.asarray(line.get_xdata()).tobytes(),np.asarray(line.get_ydata()).tobytes()])
        for col in ax.collections:
            out.append(np.asarray(col.get_offsets()).tobytes())
            if hasattr(col,'get_segments'):
                out.extend(np.asarray(x).tobytes() for x in col.get_segments())
    return out

def save(fig,out,name):
    before=geometry(fig)
    if name=='Figure_3':
        handles=[Line2D([0],[0],marker='o',lw=.8,ms=4,color=M.ORANGE,label='Source q < 0.05'),
                 Line2D([0],[0],marker='o',lw=.8,ms=4,color=M.GREY,label='Source q ≥ 0.05')]
        fig.legend(handles=handles,loc='center',bbox_to_anchor=(.515,.439),ncol=2,frameon=False,
                   fontsize=8,handlelength=1.4,columnspacing=1.4,handletextpad=.5)
    if name=='Figure_5':
        for txt in list(fig.texts):
            if txt.get_text().startswith('●'):
                txt.remove()
        handles=[Line2D([0],[0],marker='o',lw=.8,ms=3,color=M.BLUE,label='Unadjusted'),
                 Line2D([0],[0],marker='o',lw=.8,ms=3,color=M.ORANGE,label='Age and centre adjusted')]
        fig.legend(handles=handles,loc='center',bbox_to_anchor=(.595,.913),ncol=2,frameon=False,
                   fontsize=8,handlelength=1.4,columnspacing=1.4,handletextpad=.5)
        for ax in fig.axes[-2:]:
            for col in ax.collections:
                col.set_facecolor('#4D4D4D')
    assert geometry(fig)==before, 'Display changes altered plotted data geometry'
    crop={'Figure_1':1.45,'Figure_3':.375,'Figure_5':.25}.get(name,0)
    bbox=Bbox.from_extents(0,crop,*fig.get_size_inches()) if crop else None
    label='Supplementary_Figure_8' if name=='Figure_5' else name
    for ext in ['png','pdf','svg']:
        fig.savefig(out/f'{label}.{ext}',format=ext,dpi=400,bbox_inches=bbox)
    M.plt.close(fig)

def figure1(out):
    fig=M.plt.figure(figsize=(7.5,4.4))
    ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.5,.95,'Cohort consistency and cellular interpretation',ha='center',fontsize=13,weight='bold')
    cards=[
      (.025,M.BLUE,'Blood','Which immune signals recur\nacross cohorts?',
       'Acute GBS and healthy controls\nWhole blood, leukocytes, monocytes\n3 cohorts, 31 participants\n7 gene modules'),
      (.355,M.TEAL,'Cerebrospinal fluid','Do published proteins show\nrelated signals?',
       'GBS, CIDP and healthy controls\n11 paired GBS protein estimates\n16 disease comparison estimates\nSource tests and pooled uncertainty'),
      (.685,M.ORANGE,'Peripheral nerve','How does composition affect\nexpression estimates?',
       '9 CIDP and 11 CIAP donors\nMacrophage Fc and composition\nExternal nerve: 4 CIDP, 9 VN\nWhole tissue Fc comparison')]
    for x,c,h,sub,body in cards:
        ax.add_patch(FancyBboxPatch((x,.37),.29,.47,boxstyle='round,pad=0.01,rounding_size=.012',fc='#f6f8f9',ec='#d8e0e4',lw=.8))
        ax.add_patch(FancyBboxPatch((x,.74),.29,.10,boxstyle='round,pad=0.01,rounding_size=.012',fc=c,ec=c))
        ax.text(x+.145,.79,h,color='white',weight='bold',fontsize=12,ha='center',va='center')
        ax.text(x+.145,.675,sub,ha='center',va='center',weight='bold',fontsize=8.7,linespacing=1.3)
        ax.text(x+.145,.50,body,ha='center',va='center',fontsize=8.2,linespacing=1.9)
    save(fig,out,'Figure_1')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'figures')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    M.save=save
    figure1(args.output)
    M.figure3(ROOT/'data',args.output)
    M.figure5(ROOT/'data',ROOT/'data/nerve_clinical_plot_data.csv',args.output)
    print('Rendered Figure 1, Figure 3 and Supplementary Figure S8. Data geometry unchanged.')
