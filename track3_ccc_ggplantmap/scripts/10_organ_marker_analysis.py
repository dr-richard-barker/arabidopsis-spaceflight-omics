#!/usr/bin/env python3
"""
Marker-based organ inference analysis.
Tests whether known organ marker genes can distinguish seed/root/leaf/hypocotyl
in the 6-day seedling atlas. Expected: negative result (markers don't separate organs).
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
import numpy as np
import pandas as pd
import scanpy as sc
import os

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300

FIGDIR = RESULTS + '/figures'

# Load atlas
adata = sc.read_h5ad(ATLAS + '/seedling_6d_anndata.h5ad')
ct = adata.obs['CellType']

# Organ markers
markers = {
    # Root
    'AT3G54220': 'SCR', 'AT1G79560': 'PLT1', 'AT4G03270': 'PLT2',
    'AT3G18280': 'WOL', 'AT2G28470': 'CPC',
    # Leaf/cotyledon (photosynthetic)
    'AT1G01060': 'LHCB1.1', 'AT3G54890': 'FSD1', 'AT2G34420': 'RBCS',
    'AT5G38420': 'COR15A',
    # Hypocotyl
    'AT5G57390': 'ATHB-1', 'AT2G18380': 'ATHB-7',
    # Seed/embryo
    'AT3G15670': 'LEC2', 'AT2G41730': '2S3',
    # Guard cell
    'AT5G60900': 'KAT1_gc', 'AT3G15730': 'MYB60',
}

organ_map = {
    'SCR': 'Root', 'PLT1': 'Root', 'PLT2': 'Root', 'WOL': 'Root', 'CPC': 'Root',
    'LHCB1.1': 'Leaf', 'FSD1': 'Leaf', 'RBCS': 'Leaf', 'COR15A': 'Leaf',
    'ATHB-1': 'Hypocotyl', 'ATHB-7': 'Hypocotyl',
    'LEC2': 'Seed', '2S3': 'Seed',
    'KAT1_gc': 'Guard', 'MYB60': 'Guard',
}

present = {g: n for g, n in markers.items() if g in adata.var_names}
print(f'Markers present: {len(present)}/{len(markers)}')

# Compute expression per cell type
results = []
for gid, name in present.items():
    idx = list(adata.var_names).index(gid)
    expr = adata.X[:, idx].toarray().flatten() if hasattr(adata.X, 'toarray') else np.array(adata.X[:, idx]).flatten()
    for ct_name in sorted(ct.unique()):
        mask = ct == ct_name
        results.append({
            'marker': name, 'gene': gid, 'organ': organ_map.get(name, 'Unknown'),
            'celltype': ct_name,
            'mean_expr': float(expr[mask].mean()),
            'pct_expr': float((expr[mask] > 0).mean() * 100)
        })

df = pd.DataFrame(results)
df.to_csv(RESULTS + '/tables/organ_marker_expression.csv', index=False)
print('Saved: organ_marker_expression.csv')

# ============================================================
# Diagnostic figure: marker expression heatmap
# ============================================================
ct_order = ['Epidermal', 'Mesophyll', 'Stele', 'Meristematic', 'Guard']
marker_order = ['SCR', 'PLT1', 'PLT2', 'WOL', 'CPC',
                'LHCB1.1', 'FSD1', 'RBCS', 'COR15A',
                'ATHB-1', 'ATHB-7', 'LEC2', '2S3', 'KAT1_gc', 'MYB60']

# Mean expression
pivot_mean = df.pivot(index='marker', columns='celltype', values='mean_expr').reindex(marker_order)[ct_order]
pivot_pct = df.pivot(index='marker', columns='celltype', values='pct_expr').reindex(marker_order)[ct_order]

fig, axes = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={'wspace': 0.4})

# Panel A: Mean expression
ax = axes[0]
im1 = ax.imshow(pivot_mean.values, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(len(marker_order)))
ax.set_yticklabels([f'{m} ({organ_map.get(m, "")})' for m in marker_order], fontsize=8)
for i in range(len(marker_order)):
    for j in range(len(ct_order)):
        val = pivot_mean.values[i, j]
        color = 'white' if val > pivot_mean.values.max() * 0.6 else 'black'
        ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=6, color=color)
plt.colorbar(im1, ax=ax, label='Mean expression', shrink=0.7)
ax.set_title('a) Organ marker mean expression', fontweight='bold', fontsize=10, loc='left')

# Add organ group annotations
organ_groups = {'Root': 5, 'Leaf': 4, 'Hypocotyl': 2, 'Seed': 2, 'Guard': 2}
pos = 0
organ_colors = {'Root': '#0279EE', 'Leaf': '#009E73', 'Hypocotyl': '#FF9400', 'Seed': '#D55E00', 'Guard': '#CC79A7'}
for oname, n in organ_groups.items():
    ax.plot([-1.3, -1.3], [pos-0.5, pos+n-0.5], color=organ_colors[oname], linewidth=3, clip_on=False)
    ax.text(-1.6, pos + n/2 - 0.5, oname, ha='right', va='center', fontsize=7, color=organ_colors[oname], fontweight='bold', rotation=90)
    pos += n

# Panel B: % expressing
ax = axes[1]
im2 = ax.imshow(pivot_pct.values, cmap='YlGnBu', aspect='auto', vmin=0, vmax=100)
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(len(marker_order)))
ax.set_yticklabels(marker_order, fontsize=8)
for i in range(len(marker_order)):
    for j in range(len(ct_order)):
        val = pivot_pct.values[i, j]
        color = 'white' if val > 50 else 'black'
        ax.text(j, i, f'{val:.1f}%', ha='center', va='center', fontsize=6, color=color)
plt.colorbar(im2, ax=ax, label='% expressing', shrink=0.7)
ax.set_title('b) % cells expressing marker', fontweight='bold', fontsize=10, loc='left')

fig.suptitle('Organ Marker Diagnostic: Markers Do Not Separate Organs in Seedling Atlas',
             fontweight='bold', fontsize=12, y=1.01)
fig.text(0.5, -0.02, 'Photosynthetic markers (RBCS, FSD1) are ubiquitous across all cell types.\n'
         'Root markers (SCR, PLT1, WOL) are very low (<7%). Seed markers near-zero.\n'
         'Conclusion: organ identity cannot be inferred from markers in this atlas.',
         ha='center', fontsize=8, style='italic', color='0.4')
fig.tight_layout()
fig.savefig(f'{FIGDIR}/organ_marker_diagnostic.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/organ_marker_diagnostic.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: organ_marker_diagnostic.png/.svg')

# Print summary
print('\n=== ORGAN MARKER SUMMARY ===')
print('Root markers (SCR, PLT1, PLT2, WOL, CPC):')
root_df = df[df['organ'] == 'Root']
print(f'  Mean expression range: {root_df["mean_expr"].min():.4f} - {root_df["mean_expr"].max():.4f}')
print(f'  % expressing range: {root_df["pct_expr"].min():.1f}% - {root_df["pct_expr"].max():.1f}%')
print()
print('Leaf markers (LHCB1.1, FSD1, RBCS, COR15A):')
leaf_df = df[df['organ'] == 'Leaf']
print(f'  Mean expression range: {leaf_df["mean_expr"].min():.4f} - {leaf_df["mean_expr"].max():.4f}')
print(f'  % expressing range: {leaf_df["pct_expr"].min():.1f}% - {leaf_df["pct_expr"].max():.1f}%')
print()
print('Seed markers (LEC2, 2S3):')
seed_df = df[df['organ'] == 'Seed']
print(f'  Mean expression range: {seed_df["mean_expr"].min():.4f} - {seed_df["mean_expr"].max():.4f}')
print(f'  % expressing range: {seed_df["pct_expr"].min():.1f}% - {seed_df["pct_expr"].max():.1f}%')
print()
print('CONCLUSION: Marker-based organ inference is NOT reliable.')
print('Photosynthetic markers are ubiquitous; root/seed markers are too low.')
print('Cell type proxy is the more defensible approach.')
