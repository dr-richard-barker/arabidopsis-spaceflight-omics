"""
Build Figure 8 panel B: CBL9-CIPK23-AKT1 molecular cascade Sankey
Using matplotlib bezier patches for full label control (no clipping).
"""
# --- portable paths (de-sandboxed; replaces /mnt/results and /workspace) ---
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
REPO_ROOT = _os.path.abspath(_os.path.join(_HERE, '..', '..'))
RESULTS = _os.environ.get("ASO_ROOT", REPO_ROOT)          # holds tables/ and figures/
ATLAS   = _os.environ.get("ASO_ATLAS", _os.path.join(REPO_ROOT, "atlas"))  # large intermediates (not shipped)
WORK    = _os.environ.get("ASO_WORK", _os.path.join(REPO_ROOT, "work"))    # scratch outputs
_os.makedirs(WORK, exist_ok=True)
# --- end portable paths ---

import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

# ---- Data ----
expr = pd.read_csv(ATLAS + '/ca2_k_expr_per_celltype.csv')
coexpr_guard = pd.read_csv(RESULTS + '/tables/ca2_k_coexpression_Guard.csv', index_col=0)
coexpr_guard.index = [x.replace('_Guard', '') for x in coexpr_guard.index]

def get_guard_expr(alias):
    sub = expr[(expr['alias'] == alias) & (expr['celltype'] == 'Guard')]
    return sub.iloc[0]['mean_expr'] if len(sub) > 0 else 0.01

def get_coexpr(g1, g2):
    if g1 in coexpr_guard.index and g2 in coexpr_guard.columns:
        return coexpr_guard.loc[g1, g2]
    return 0.0

# ---- Node definitions (column, row position) ----
# Columns: 0=Ca2+, 1=sensors, 2=kinases, 3=channels/targets, 4=output
nodes = {
    'Ca\u00b2\u207a':  {'col': 0, 'row': 2.5, 'color': '#E9ED4C',  'expr': 1.0,   'type': 'ligand'},
    'CBL1':   {'col': 1, 'row': 4.0, 'color': '#0279EE',  'expr': get_guard_expr('CBL1'),   'type': 'sensor'},
    'CBL9':   {'col': 1, 'row': 3.0, 'color': '#0279EE',  'expr': get_guard_expr('CBL9'),   'type': 'sensor'},
    'CaM3':   {'col': 1, 'row': 2.0, 'color': '#FD9BED',  'expr': get_guard_expr('CaM3'),   'type': 'sensor'},
    'CML24':  {'col': 1, 'row': 1.0, 'color': '#FD9BED',  'expr': get_guard_expr('CML24'),  'type': 'sensor'},
    'CIPK23': {'col': 2, 'row': 3.5, 'color': '#75A025',  'expr': get_guard_expr('CIPK23'), 'type': 'kinase'},
    'CDPK4':  {'col': 2, 'row': 2.0, 'color': '#FF9400',  'expr': get_guard_expr('CDPK4'),  'type': 'kinase'},
    'CDPK3':  {'col': 2, 'row': 1.0, 'color': '#FF9400',  'expr': get_guard_expr('CDPK3'),  'type': 'kinase'},
    'AKT1':   {'col': 3, 'row': 3.5, 'color': '#FF9400',  'expr': get_guard_expr('AKT1'),   'type': 'channel'},
    'KC1':    {'col': 3, 'row': 1.5, 'color': '#009E73',  'expr': get_guard_expr('KC1'),    'type': 'channel'},
    'ANNAT1': {'col': 3, 'row': 2.5, 'color': '#56B4E9',  'expr': get_guard_expr('ANNAT1'), 'type': 'target'},
    'K\u207a uptake': {'col': 4, 'row': 2.5, 'color': '#ECE9E2', 'expr': 0.3, 'type': 'output'},
}

# ---- Edge definitions ----
edges = [
    # Ca2+ -> sensors
    ('Ca\u00b2\u207a', 'CBL1',   '#0279EE', 0.4),
    ('Ca\u00b2\u207a', 'CBL9',   '#0279EE', 0.4),
    ('Ca\u00b2\u207a', 'CaM3',   '#FD9BED', 0.4),
    ('Ca\u00b2\u207a', 'CML24',  '#FD9BED', 0.4),
    # sensors -> kinases
    ('CBL1',   'CIPK23', '#0279EE', 0.4),
    ('CBL9',   'CIPK23', '#0279EE', 0.4),
    ('CaM3',   'CDPK4',  '#FD9BED', 0.4),
    ('CaM3',   'CDPK3',  '#FD9BED', 0.4),
    # kinases -> channels
    ('CIPK23', 'AKT1',   '#FF9400', 0.7),  # crosstalk - thicker, orange
    ('CDPK4',  'ANNAT1', '#56B4E9', 0.65),
    ('CML24',  'KC1',    '#FD9BED', 0.4),
    # channels -> output
    ('AKT1',   'K\u207a uptake', '#FF9400', 0.6),
    ('KC1',    'K\u207a uptake', '#009E73', 0.4),
]

# ---- Layout ----
COL_X = {0: 0.08, 1: 0.28, 2: 0.52, 3: 0.72, 4: 0.92}
ROW_SCALE = 0.14  # row units to figure fraction

def node_xy(name):
    n = nodes[name]
    return COL_X[n['col']], n['row'] * ROW_SCALE

def node_height(name):
    """Node height proportional to expression"""
    return max(nodes[name]['expr'] * 0.08, 0.012)

# ---- Draw ----
fig, ax = plt.subplots(figsize=(10, 9))
ax.set_xlim(0, 1)
ax.set_ylim(0, 0.68)
ax.axis('off')

# Draw bezier flow bands
from matplotlib.path import Path
import matplotlib.patches as patches

def draw_flow(ax, x0, y0_center, x1, y1_center, width, color, alpha=0.45):
    """Draw a bezier flow band between two nodes."""
    hw = width / 2
    # Control points for smooth S-curve
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0_center - hw),
        (cx, y0_center - hw),
        (cx, y1_center - hw),
        (x1, y1_center - hw),
        (x1, y1_center + hw),
        (cx, y1_center + hw),
        (cx, y0_center + hw),
        (x0, y0_center + hw),
        (x0, y0_center - hw),
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, edgecolor='none', alpha=alpha, zorder=1)
    ax.add_patch(patch)

# Draw flows
for (src, tgt, color, alpha) in edges:
    x0, y0 = node_xy(src)
    x1, y1 = node_xy(tgt)
    # Flow width proportional to target expression
    w = max(nodes[tgt]['expr'] * 0.025, 0.004)
    draw_flow(ax, x0, y0, x1, y1, w, color, alpha)

# Draw nodes (rectangles)
NODE_W = 0.025
for name, nd in nodes.items():
    x, y = node_xy(name)
    h = node_height(name)
    rect = mpatches.FancyBboxPatch(
        (x - NODE_W/2, y - h/2), NODE_W, h,
        boxstyle='round,pad=0.002',
        facecolor=nd['color'], edgecolor='#333333', linewidth=1.2, zorder=3
    )
    ax.add_patch(rect)

# Draw node labels
LABEL_FONTSIZE = 11
for name, nd in nodes.items():
    x, y = node_xy(name)
    col = nd['col']
    # Left-side nodes: label to the left; right-side: label to the right; middle: centered above
    if col == 0:
        ax.text(x - NODE_W/2 - 0.01, y, name, ha='right', va='center',
                fontsize=LABEL_FONTSIZE, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.8, edgecolor='none'))
    elif col == 4:
        ax.text(x + NODE_W/2 + 0.01, y, name, ha='left', va='center',
                fontsize=LABEL_FONTSIZE, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.8, edgecolor='none'))
    else:
        # Middle columns: label above node
        ax.text(x, y + node_height(name)/2 + 0.012, name, ha='center', va='bottom',
                fontsize=LABEL_FONTSIZE - 1, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.08', facecolor='white', alpha=0.85, edgecolor='none'))

# ---- Co-expression r-value annotations on key edges ----
key_coexpr = [
    ('CDPK4', 'ANNAT1'),
    ('CIPK23', 'AKT1'),
    ('CaM3', 'CDPK4'),
]
for g1, g2 in key_coexpr:
    r = get_coexpr(g1, g2)
    if abs(r) >= 0.1:
        x0, y0 = node_xy(g1)
        x1, y1 = node_xy(g2)
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        ax.text(xm, ym + 0.015, f'r={r:.2f}', ha='center', va='bottom',
                fontsize=9.5, fontstyle='italic', fontweight='bold', color='#222222',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='#FFFDE7', alpha=0.95, edgecolor='#CCAA00', linewidth=0.8))

# ---- Column headers ----
col_labels = {0: 'Signal', 1: 'Sensors', 2: 'Kinases', 3: 'Channels', 4: 'Output'}
col_colors = {0: '#888888', 1: '#0279EE', 2: '#FF9400', 3: '#009E73', 4: '#888888'}
for col, label in col_labels.items():
    ax.text(COL_X[col], 0.66, label, ha='center', va='top',
            fontsize=10, fontweight='bold', color=col_colors[col],
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.9,
                      edgecolor=col_colors[col], linewidth=1))

# ---- Legend ----
legend_items = [
    mpatches.Patch(color='#E9ED4C', label='Ca\u00b2\u207a signal'),
    mpatches.Patch(color='#0279EE', label='CBL sensor'),
    mpatches.Patch(color='#FD9BED', label='CaM/CML sensor'),
    mpatches.Patch(color='#75A025', label='CIPK kinase'),
    mpatches.Patch(color='#FF9400', label='CDPK kinase'),
    mpatches.Patch(color='#009E73', label='K\u207a channel'),
    mpatches.Patch(color='#56B4E9', label='Ca\u00b2\u207a target'),
    mpatches.Patch(color='#ECE9E2', label='K\u207a uptake', edgecolor='grey'),
]
ax.legend(handles=legend_items, loc='lower center', ncol=4,
          fontsize=8.5, framealpha=0.9, edgecolor='grey',
          bbox_to_anchor=(0.5, -0.02))

# ---- Title ----
ax.set_title('CBL9\u2013CIPK23\u2013AKT1 molecular cascade\n'
             'Guard cell expression (node height) \u00b7 co-expression r (italic labels)',
             fontsize=13, fontweight='bold', pad=12)

plt.tight_layout(pad=0.5)
out = WORK + '/figure8b_matplotlib.png'
plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

from PIL import Image
import os
img = Image.open(out)
print(f'Saved: {out} - {img.size} - {os.path.getsize(out):,} bytes')
