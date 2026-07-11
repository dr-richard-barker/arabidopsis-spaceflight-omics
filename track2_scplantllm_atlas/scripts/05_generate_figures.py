#!/usr/bin/env python3
"""
Figure 4: LASSO gene expression on plant anatomy (from R composites)
Figure 5: Ca2+ signaling pathway (ggpathway + conceptual schematic)
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

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300

FIGDIR = RESULTS + '/figures'

# ============================================================
# FIGURE 4: LASSO genes on leaf + root (from R composites)
# ============================================================
print('=== Figure 4: LASSO genes on plant anatomy ===')

fig, axes = plt.subplots(2, 1, figsize=(12, 14))

# Panel a: Leaf composite
leaf_img = mpimg.imread(f'{FIGDIR}/ggplantmap_composite1_lasso_leaf.png')
axes[0].imshow(leaf_img)
axes[0].set_title('a) LASSO biomarker gene expression on leaf cross-section',
                  fontweight='bold', loc='left', fontsize=11)
axes[0].axis('off')

# Panel b: Root composite
root_img = mpimg.imread(f'{FIGDIR}/ggplantmap_composite2_lasso_root.png')
axes[1].imshow(root_img)
axes[1].set_title('b) LASSO biomarker gene expression on root tip cross-section',
                  fontweight='bold', loc='left', fontsize=11)
axes[1].axis('off')

fig.suptitle('Spatial Expression of LASSO Biomarker Genes on Arabidopsis Anatomy',
             fontsize=13, fontweight='bold', y=0.98)
fig.tight_layout()
fig.savefig(f'{FIGDIR}/Figure4_ggPlantmap_LASSO_genes.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/Figure4_ggPlantmap_LASSO_genes.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print('  Saved: Figure4_ggPlantmap_LASSO_genes.png/.svg')

# ============================================================
# FIGURE 5: Ca2+ pathway (ggpathway + conceptual schematic)
# ============================================================
print('=== Figure 5: Ca2+ signaling pathway ===')

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Panel a: ggpathway network
ggpath_img = mpimg.imread(f'{FIGDIR}/ccc_ca2_pathway_ggpathway.png')
axes[0].imshow(ggpath_img)
axes[0].set_title('a) Ca2+ signaling network (data-driven, ggpathway)',
                  fontweight='bold', loc='left', fontsize=11)
axes[0].axis('off')

# Panel b: Conceptual schematic
schematic_img = mpimg.imread(f'{FIGDIR}/ccc_ca2_conceptual_schematic.png')
axes[1].imshow(schematic_img)
axes[1].set_title('b) Ca2+ signaling cascade under microgravity (conceptual)',
                  fontweight='bold', loc='left', fontsize=11)
axes[1].axis('off')

fig.suptitle('Ca2+ Signaling Pathway: Dominant Cell-Cell Communication Axis',
             fontsize=13, fontweight='bold', y=0.98)
fig.tight_layout()
fig.savefig(f'{FIGDIR}/Figure5_Ca2_pathway.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/Figure5_Ca2_pathway.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print('  Saved: Figure5_Ca2_pathway.png/.svg')

print('\n=== All 6 manuscript figures complete ===')
for i in range(1, 7):
    names = {
        1: 'Figure1_LASSO_panel',
        2: 'Figure2_scPlantLLM_clustering',
        3: 'Figure3_CCC_communication',
        4: 'Figure4_ggPlantmap_LASSO_genes',
        5: 'Figure5_Ca2_pathway',
        6: 'Figure6_integration'
    }
    name = names[i]
    png_path = f'{FIGDIR}/{name}.png'
    if os.path.exists(png_path):
        size = os.path.getsize(png_path) / 1024
        print(f'  Figure {i}: {name}.png ({size:.0f} KB)')
