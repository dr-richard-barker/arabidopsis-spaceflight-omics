"""
Step 3: Reassemble Figure 5 with the refined 20-node Ca2+ pathway network.
Panel a: refined ccc_ca2_pathway_ggpathway.png (20 nodes, stress-majorization)
Panel b: ccc_ca2_conceptual_schematic.png (existing conceptual schematic)
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

fig, axes = plt.subplots(1, 2, figsize=(20, 9), gridspec_kw={'width_ratios': [1, 1.2]})
fig.patch.set_facecolor('white')

# Panel a: refined 20-node network
ggpath_img = mpimg.imread(f'{fig_dir}/ccc_ca2_pathway_ggpathway.png')
axes[0].imshow(ggpath_img)
axes[0].set_axis_off()
axes[0].set_title('a) Data-driven Ca2+ signaling network (20 nodes, stress-majorization layout)',
                  fontsize=13, fontweight='bold', pad=10)

# Panel b: conceptual schematic
schematic_img = mpimg.imread(f'{fig_dir}/ccc_ca2_conceptual_schematic.png')
axes[1].imshow(schematic_img)
axes[1].set_axis_off()
axes[1].set_title('b) Ca2+ signaling cascade under microgravity (conceptual schematic)',
                  fontsize=13, fontweight='bold', pad=10)

plt.tight_layout()

out_png = f'{fig_dir}/Figure5_Ca2_pathway.png'
out_svg = f'{fig_dir}/Figure5_Ca2_pathway.svg'
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
plt.close()

# Verify
img = Image.open(out_png)
arr = np.array(img)
nw = (arr < 250).any(axis=2).sum()
print(f"Figure 5 (reassembled): {img.size}, {os.path.getsize(out_png)/1024:.0f} KB, non-white={100*nw/(arr.shape[0]*arr.shape[1]):.1f}%")
print(f"SVG: {os.path.getsize(out_svg)/1024:.0f} KB")
