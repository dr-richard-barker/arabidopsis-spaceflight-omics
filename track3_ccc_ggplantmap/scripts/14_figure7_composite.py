"""
Step 7: Composite Figure 7 — Ca2+/K+ signaling circuit.
4 panels: (A) Ca2+/K+ pathway network, (B) co-expression heatmap,
(C) ggPlantmap seedling spatial map, (D) CCC circuit diagram.
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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch
import os

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'

fig_dir = RESULTS + '/figures'

panels = {
    'A': ('ccc_ca2_k_pathway_ggpathway.png', 'Ca2+/K+ pathway network'),
    'B': ('ca2_k_coexpression_heatmap.png', 'Co-expression by cell type'),
    'C': ('ggplantmap_ca2_k_circuit_seedling.png', 'Spatial expression (seedling map)'),
    'D': ('ca2_k_ccc_circuit.png', 'Cell-to-cell CCC circuit'),
}

# Layout: 2x2 grid
fig, axes = plt.subplots(2, 2, figsize=(18, 18))
fig.patch.set_facecolor('white')

panel_labels = {
    (0,0): 'A', (0,1): 'B',
    (1,0): 'C', (1,1): 'D',
}

for (r, c), label in panel_labels.items():
    ax = axes[r][c]
    fname, title = panels[label]
    img = mpimg.imread(os.path.join(fig_dir, fname))
    ax.imshow(img)
    ax.set_axis_off()
    # Panel label
    ax.text(-0.02, 1.02, label, transform=ax.transAxes,
            fontsize=20, fontweight='bold', va='bottom', ha='right')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=10)

fig.suptitle('Ca2+/K+ signaling circuit in Arabidopsis seedlings under spaceflight',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
out_png = RESULTS + '/figures/Figure7_ca2_k_circuit_composite.png'
out_svg = RESULTS + '/figures/Figure7_ca2_k_circuit_composite.svg'
fig.savefig(out_png, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
plt.close()

# Verify
from PIL import Image
import numpy as np
img = Image.open(out_png)
arr = np.array(img)
nw = (arr < 250).any(axis=2).sum()
print(f"Figure7: {img.size}, {os.path.getsize(out_png)/1024:.0f} KB, non-white={100*nw/(arr.shape[0]*arr.shape[1]):.1f}%")
