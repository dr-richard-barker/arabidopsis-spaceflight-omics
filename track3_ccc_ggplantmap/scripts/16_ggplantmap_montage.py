"""
Step 2: ggPlantmap montage — Supplementary Figure.
Vertical 3-panel stack of seedling, root tip, and leaf ggPlantmap Ca2+/K+ circuit maps.
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
import os
import numpy as np
from PIL import Image

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'

fig_dir = RESULTS + '/figures'

panels = [
    ('ggplantmap_ca2_k_circuit_seedling.png', 'a) Seedling (cotyledon, hypocotyl, root)',
     'Whole-plant expression of 6 circuit genes: AKT1, CBL9, CIPK23, CML24, KC1, CDPK3'),
    ('ggplantmap_ca2_k_circuit_roottip.png', 'b) Root tip longitudinal section',
     'Ca2+ cascade gene expression across 12 root zones (apical meristem → elongation → maturation)'),
    ('ggplantmap_ca2_k_circuit_leaf.png', 'c) Leaf cross-section',
     'K+ channel distribution across 8 leaf regions (epidermis, palisade, spongy mesophyll, vascular)'),
]

fig, axes = plt.subplots(3, 1, figsize=(14, 20))
fig.patch.set_facecolor('white')

for idx, (fname, title, subtitle) in enumerate(panels):
    ax = axes[idx]
    img = mpimg.imread(os.path.join(fig_dir, fname))
    ax.imshow(img)
    ax.set_axis_off()

    # Panel label
    ax.text(-0.01, 1.02, title, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='bottom', ha='left')
    ax.text(-0.01, 0.98, subtitle, transform=ax.transAxes,
            fontsize=10, style='italic', va='top', ha='left', color='#555555')

fig.suptitle('Spatial expression of Ca2+/K+ circuit genes across Arabidopsis anatomy',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])

out_png = RESULTS + '/figures/SupplementaryFigure_ggplantmap_montage.png'
out_svg = RESULTS + '/figures/SupplementaryFigure_ggplantmap_montage.svg'
fig.savefig(out_png, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
plt.close()

# Verify
img = Image.open(out_png)
arr = np.array(img)
nw = (arr < 250).any(axis=2).sum()
print(f"Supp Fig (ggPlantmap montage): {img.size}, {os.path.getsize(out_png)/1024:.0f} KB, non-white={100*nw/(arr.shape[0]*arr.shape[1]):.1f}%")
print(f"SVG: {os.path.getsize(out_svg)/1024:.0f} KB")
