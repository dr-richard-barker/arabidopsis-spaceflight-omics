#!/usr/bin/env python3
"""
Co-expression analysis of Ca2+/K+ pathway genes.
Creates: (1) 25x25 co-expression heatmap, (2) per-cell-type expression heatmap.
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
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import os

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300

FIGDIR = RESULTS + '/figures'

# Load data
corr = pd.read_csv(RESULTS + '/tables/ca2_k_coexpression_matrix.csv', index_col=0)
expr = pd.read_csv(ATLAS + '/ca2_k_expr_per_celltype.csv')

# ============================================================
# Figure 1: Co-expression heatmap (25x25)
# ============================================================
print('=== Co-expression heatmap ===')

# Order genes by type for visual grouping
gene_order = ['CBL1','CBL9','CBL2','CIPK23','CIPK1','CIPK9',
              'CML38','CML42','CML24','CaM7','CaM3',
              'CDPK3','CDPK4','CAX3','NRT2.1','SOS1','ANNAT1',
              'AKT1','AKT2','GORK','KAT1','KAT2','KC1',
              'HSP22.0','TRAF-like']
corr_ordered = corr.loc[gene_order, gene_order]

fig, ax = plt.subplots(figsize=(12, 10))
# Diverging colormap centered at 0
cmap = plt.cm.RdBu_r
vmax = max(abs(corr_ordered.values.min()), abs(corr_ordered.values.max()))
im = ax.imshow(corr_ordered.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect='auto')

ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_order, rotation=45, ha='right', fontsize=8)
ax.set_yticks(range(len(gene_order)))
ax.set_yticklabels(gene_order, fontsize=8)

# Add correlation values
for i in range(len(gene_order)):
    for j in range(len(gene_order)):
        val = corr_ordered.values[i, j]
        color = 'white' if abs(val) > vmax * 0.6 else 'black'
        ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=5.5, color=color)

# Add type group annotations on top
type_groups = {
    'CBL sensors': ['CBL1','CBL9','CBL2'],
    'CIPK kinases': ['CIPK23','CIPK1','CIPK9'],
    'CaM/CML': ['CML38','CML42','CML24','CaM7','CaM3'],
    'CDPKs': ['CDPK3','CDPK4'],
    'Ca2+ targets': ['CAX3','NRT2.1','SOS1','ANNAT1'],
    'K+ channels': ['AKT1','AKT2','GORK','KAT1','KAT2','KC1'],
    'LASSO': ['HSP22.0','TRAF-like'],
}
type_colors = {
    'CBL sensors': '#0279EE', 'CIPK kinases': '#0279EE',
    'CaM/CML': '#FD9BED', 'CDPKs': '#FF9400',
    'Ca2+ targets': '#56B4E9', 'K+ channels': '#009E73',
    'LASSO': '#D55E00',
}
pos = 0
for gname, genes in type_groups.items():
    n = len(genes)
    ax.plot([pos-0.5, pos+n-0.5], [-1.2, -1.2], color=type_colors[gname], linewidth=3, clip_on=False)
    ax.text(pos + n/2 - 0.5, -1.8, gname, ha='center', va='bottom', fontsize=7, color=type_colors[gname], fontweight='bold')
    pos += n

plt.colorbar(im, ax=ax, label='Pearson r', shrink=0.7)
ax.set_title('Ca2+/K+ Pathway Gene Co-expression (all 9,998 cells)', fontweight='bold', fontsize=12, pad=20)
fig.tight_layout()
fig.savefig(f'{FIGDIR}/ca2_k_coexpression_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/ca2_k_coexpression_heatmap.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: ca2_k_coexpression_heatmap.png/.svg')

# ============================================================
# Figure 2: Per-cell-type expression heatmap
# ============================================================
print('=== Per-cell-type expression heatmap ===')

# Pivot: genes x cell types (mean expression)
expr_ct = expr[expr['celltype'] != 'ALL'].pivot(index='alias', columns='celltype', values='mean_expr')
# Reorder cell types
ct_order = ['Epidermal', 'Mesophyll', 'Stele', 'Meristematic', 'Guard']
expr_ct = expr_ct[ct_order]
# Reorder genes by type
expr_ct = expr_ct.loc[gene_order]

fig, axes = plt.subplots(1, 2, figsize=(14, 10), gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.4})

# Panel A: Mean expression
ax = axes[0]
im1 = ax.imshow(expr_ct.values, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(len(gene_order)))
ax.set_yticklabels(gene_order, fontsize=8)
for i in range(len(gene_order)):
    for j in range(len(ct_order)):
        val = expr_ct.values[i, j]
        color = 'white' if val > expr_ct.values.max() * 0.6 else 'black'
        ax.text(j, i, f'{val:.3f}', ha='center', va='center', fontsize=6, color=color)
plt.colorbar(im1, ax=ax, label='Mean expression', shrink=0.7)
ax.set_title('a) Mean expression per cell type', fontweight='bold', fontsize=10, loc='left')

# Panel B: % expressing
expr_pct = expr[expr['celltype'] != 'ALL'].pivot(index='alias', columns='celltype', values='pct_expr')
expr_pct = expr_pct[ct_order].loc[gene_order]
ax = axes[1]
im2 = ax.imshow(expr_pct.values, cmap='YlGnBu', aspect='auto', vmin=0, vmax=100)
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(len(gene_order)))
ax.set_yticklabels(gene_order, fontsize=8)
for i in range(len(gene_order)):
    for j in range(len(ct_order)):
        val = expr_pct.values[i, j]
        color = 'white' if val > 50 else 'black'
        ax.text(j, i, f'{val:.1f}%', ha='center', va='center', fontsize=6, color=color)
plt.colorbar(im2, ax=ax, label='% expressing', shrink=0.7)
ax.set_title('b) % cells expressing per cell type', fontweight='bold', fontsize=10, loc='left')

fig.suptitle('Ca2+/K+ Pathway Gene Expression Across Cell Types', fontweight='bold', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(f'{FIGDIR}/ca2_k_expression_per_celltype.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/ca2_k_expression_per_celltype.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: ca2_k_expression_per_celltype.png/.svg')

# ============================================================
# Print top co-expression pairs per cell type
# ============================================================
print('\n=== Top co-expression pairs per cell type ===')
for ct_name in ct_order:
    ct_file = RESULTS + f'/tables/ca2_k_coexpression_{ct_name}.csv'
    if os.path.exists(ct_file):
        ct_corr = pd.read_csv(ct_file, index_col=0)
        # Extract upper triangle pairs
        pairs = []
        for i in range(len(gene_order)):
            for j in range(i+1, len(gene_order)):
                g1, g2 = gene_order[i], gene_order[j]
                if g1 in ct_corr.index and g2 in ct_corr.columns:
                    r = ct_corr.loc[g1, g2] if not isinstance(ct_corr.loc[g1, g2], pd.Series) else ct_corr.loc[g1, g2].iloc[0]
                    if not np.isnan(r):
                        pairs.append((g1, g2, r))
        pairs.sort(key=lambda x: abs(x[2]), reverse=True)
        print(f'\n{ct_name}:')
        for a, b, r in pairs[:5]:
            print(f'  {a} ~ {b}: r={r:.3f}')

print('\nDone.')
