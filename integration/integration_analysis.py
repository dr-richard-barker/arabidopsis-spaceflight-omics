#!/usr/bin/env python3
"""
Figure 6: Integration summary - combines LASSO panel score map + gene function table
"""
# --- portable paths (de-sandboxed; replaces /mnt/results and /workspace) ---
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
REPO_ROOT = _os.path.abspath(_os.path.join(_HERE, '..'))
RESULTS = _os.environ.get("ASO_ROOT", REPO_ROOT)          # holds tables/ and figures/
ATLAS   = _os.environ.get("ASO_ATLAS", _os.path.join(REPO_ROOT, "atlas"))  # large intermediates (not shipped)
WORK    = _os.environ.get("ASO_WORK", _os.path.join(REPO_ROOT, "work"))    # scratch outputs
_os.makedirs(WORK, exist_ok=True)
# --- end portable paths ---

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import os

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300

FIGDIR = RESULTS + '/figures'

# Load gene function annotations
gene_ann = pd.read_csv(RESULTS + '/tables/gene_function_annotations.csv')

fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

# Panel A: LASSO panel score on leaf (load the R-generated figure)
ax_a = fig.add_subplot(gs[0, 0])
score_img_path = f'{FIGDIR}/ggplantmap_composite5_lasso_panel_score.png'
if os.path.exists(score_img_path):
    img = mpimg.imread(score_img_path)
    # The composite has both leaf and root; crop to left half for leaf
    h, w = img.shape[:2]
    ax_a.imshow(img[:, :w//2, :])
ax_a.set_title('a) LASSO panel weighted score (leaf)', fontweight='bold', loc='left', fontsize=10)
ax_a.axis('off')

# Panel B: Ca2+ signaling on leaf (load from R composite)
ax_b = fig.add_subplot(gs[0, 1])
ca2_img_path = f'{FIGDIR}/ggplantmap_composite4_ca2_signaling.png'
if os.path.exists(ca2_img_path):
    img = mpimg.imread(ca2_img_path)
    h, w = img.shape[:2]
    ax_b.imshow(img[:, :w//2, :])
ax_b.set_title('b) Ca2+ signaling strength (leaf)', fontweight='bold', loc='left', fontsize=10)
ax_b.axis('off')

# Panel C: Gene function summary table
ax_c = fig.add_subplot(gs[1, :])
ax_c.axis('off')

# Build table data
table_data = []
for _, row in gene_ann.iterrows():
    direction = 'Flight (+)' if '+' in row['direction'] else 'Ground (-)'
    func_short = row['protein_name'][:45] + '...' if len(row['protein_name']) > 45 else row['protein_name']
    table_data.append([
        row['gene_id'],
        row['alias'][:15],
        row['coefficient'],
        direction,
        func_short
    ])

col_labels = ['Gene ID', 'Alias', 'Coefficient', 'Direction', 'Protein / Function']
table = ax_c.table(cellText=table_data, colLabels=col_labels,
                   loc='center', cellLoc='left',
                   colWidths=[0.12, 0.15, 0.1, 0.12, 0.51])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.8)

# Style header
for j in range(len(col_labels)):
    table[0, j].set_facecolor('#0279EE')
    table[0, j].set_text_props(color='white', fontweight='bold')

# Style rows
for i in range(1, len(table_data) + 1):
    for j in range(len(col_labels)):
        if i % 2 == 0:
            table[i, j].set_facecolor('#ECE9E2')

ax_c.set_title('c) Top 4 LASSO biomarker gene functions', fontweight='bold', loc='left', fontsize=10, pad=15)

fig.suptitle('Integration: LASSO Panel, Ca2+ Signaling, and Gene Functions',
             fontsize=13, fontweight='bold', y=0.98)

fig.savefig(f'{FIGDIR}/Figure6_integration.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/Figure6_integration.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: Figure6_integration.png/.svg')
