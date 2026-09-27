"""Final display changes only. Run from this directory; aggregate inputs are frozen."""
from pathlib import Path
import argparse
import importlib.util
from io import BytesIO
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('base', ROOT / 'base_make_figures.py')
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
M.plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.serif': ['DejaVu Serif']})

def geometry(fig):
    values = []
    for ax in fig.axes:
        values.extend([repr(ax.get_xlim()), repr(ax.get_ylim())])
        for line in ax.lines:
            values.extend([np.asarray(line.get_xdata()).tobytes(), np.asarray(line.get_ydata()).tobytes()])
        for col in ax.collections:
            values.append(np.asarray(col.get_offsets()).tobytes())
            if hasattr(col, 'get_segments'):
                values.extend(np.asarray(v).tobytes() for v in col.get_segments())
        for im in ax.images:
            values.append(np.asarray(im.get_array()).tobytes())
    return values

def save(fig, out, name):
    before = geometry(fig)
    if name == 'Figure_4':
        # Base Figure 4 is Figure 5 in the current manuscript.
        fig.legend(handles=[Line2D([], [], marker='D', linestyle='none', markersize=3,
                   color='#333333', label='q < 0.05 within cell group')],
                   loc='center', bbox_to_anchor=(.795,.438), frameon=False, fontsize=7.8)
        fig.legend(handles=[Line2D([], [], marker='o', linewidth=.8, markersize=3,
                   color=M.BLUE, label='Without composition'),
                   Line2D([], [], marker='o', linewidth=.8, markersize=3,
                   color=M.ORANGE, label='With composition')],
                   loc='center', bbox_to_anchor=(.755,.035), ncol=2, frameon=False,
                   fontsize=7.5, handlelength=1, columnspacing=1)
        name = 'Figure_5'
    elif name == 'Supplementary_Figure_5':
        # Replace the grid in panel B with an actual three-line table.
        ax = next(ax for ax in fig.axes if ax.tables)
        old = list(ax.tables)[0]
        rows = [[old[(i,j)].get_text().get_text() for j in range(2)] for i in range(4)]
        old.remove()
        ax.set_position([.045,.085,.43,.30])
        table = ax.table(cellText=rows, colLabels=['Gene','Published estimate'],
                         colWidths=[.24,.76], cellLoc='left', colLoc='left', bbox=[0,0,1,1])
        table.auto_set_font_size(False)
        table.set_fontsize(8.1)
        for (i,j), cell in table.get_celld().items():
            cell.set_facecolor('white')
            cell.set_edgecolor('#303030')
            cell.set_linewidth(.65)
            cell.visible_edges = 'TB' if i == 0 else ('B' if i == 4 else '')
            if i == 0:
                cell.set_text_props(weight='bold')
                cell.set_height(.12)
            else:
                cell.set_height(.22)
    elif name in ['Supplementary_Figure_6','Supplementary_Figure_7']:
        for text in list(fig.texts):
            if text.get_text().startswith(('All 15','Two cultured')):
                text.remove()
    assert geometry(fig) == before, 'A display edit changed plotted data'
    for ext in ['png','pdf','svg']:
        buffer = BytesIO()
        fig.savefig(buffer, format=ext, dpi=400, bbox_inches='tight', pad_inches=.10)
        target = out / f'{name}.{ext}'
        temporary = target.with_suffix(target.suffix + '.tmp')
        temporary.write_bytes(buffer.getvalue())
        temporary.replace(target)
    M.plt.close(fig)

def figure1(out):
    # This workflow uses reported study counts; no inferred data are plotted.
    fig = M.plt.figure(figsize=(8.4,5.1))
    ax = fig.add_axes([0,0,1,1])
    ax.set(xlim=(0,1), ylim=(0,1))
    ax.axis('off')
    columns = [
        (.018,M.BLUE,'Blood','Do recruitment signals recur\nin independent patients?',[
            ('Initial comparison','3 GBS cohorts\nWhole blood, leukocytes and monocytes'),
            ('Independent monocytes','4 GBS and 3 controls\nSame donors also assessed as PBMCs\nand sorted cell preparations'),
            ('Separate PBMC reanalysis','12 acute GBS and 5 controls\n8 paired follow up samples\nPossible cohort overlap assessed')]),
        (.350,M.TEAL,'Cerebrospinal fluid','Do published proteins show\nrelated immune changes?',[
            ('Independent source cohorts','2 GBS protein cohorts\n11 proteins with paired estimates'),
            ('Protein synthesis','Source estimates and pooled effects\nUncertainty assessed across models'),
            ('Direct disease comparisons','16 reported estimates\nGBS, CIDP and healthy controls\nWithin the replication cohort')]),
        (.682,M.ORANGE,'Peripheral nerve','How does macrophage composition\nrelate to the Fc score difference?',[
            ('Atlas macrophages','9 CIDP and 11 CIAP donors\nState composition and Fc expression'),
            ('External whole nerve','4 CIDP and 9 VN biopsies\nInternal atlas check\n9 CIDP and 5 VN donors'),
            ('Macrophage marker index','20 genes assessed in 37 atlas donors\nDonor and centre exclusions\nExternal composition adjustment')])]
    for x, colour, heading, question, boxes in columns:
        ax.add_patch(FancyBboxPatch((x,.86),.30,.095, boxstyle='round,pad=.006,rounding_size=.008',facecolor=colour,edgecolor=colour))
        ax.text(x+.15,.906,heading,color='white',weight='bold',fontsize=12,ha='center',va='center')
        ax.text(x+.15,.795,question,ha='center',va='center',weight='bold',fontsize=9.1,linespacing=1.4)
        for i,(label,body) in enumerate(boxes):
            y = .555 - .235*i
            ax.add_patch(FancyBboxPatch((x,y),.30,.19, boxstyle='round,pad=.006,rounding_size=.008',facecolor='#F7F8FA',edgecolor='#CCD4DB',linewidth=.75))
            ax.text(x+.15,y+.157,label,ha='center',va='center',fontsize=9.3,weight='bold',color=colour)
            ax.text(x+.15,y+.075,body,ha='center',va='center',fontsize=8.6,linespacing=1.5)
            if i < 2:
                ax.annotate('', xy=(x+.15,y-.036), xytext=(x+.15,y-.008),arrowprops=dict(arrowstyle='-|>',color=colour,lw=.9))
    save(fig,out,'Figure_1')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT.parents[1]/'Figures')
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    M.save = save
    figure1(args.output)
    M.figure4(ROOT/'inputs',args.output)
    M.supplement5(ROOT/'inputs',args.output)
    M.supplement6(ROOT/'inputs',args.output)
    M.supplement7(ROOT/'inputs',args.output)
    (args.output/'Figure_4D_model_rows.csv').unlink(missing_ok=True)
    print('Updated Figure 1, Figure 5, and Figures S5–S7; numerical plot geometry retained.')
