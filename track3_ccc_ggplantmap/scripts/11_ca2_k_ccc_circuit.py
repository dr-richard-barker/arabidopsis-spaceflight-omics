#!/usr/bin/env python3
"""
Cell type proxy organ mapping + CCC circuit detection.
Detects the Ca2+/K+ multicellular signaling circuit between cell types
and creates a CCC circuit diagram.
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
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300

FIGDIR = RESULTS + '/figures'

# Load CCC data
fc = pd.read_csv(RESULTS + '/tables/ccc_full_communication_table.csv')
ca2_comm = fc[fc['Signal'] == 'Ca2+']
k_comm = fc[fc['Signal'] == 'K+']

# Cell type -> organ proxy mapping
ct_to_organ = {
    'Epidermal': 'All organs (surface)',
    'Mesophyll': 'Cotyledon/leaf',
    'Stele': 'Root/hypocotyl (vascular)',
    'Meristematic': 'Root tip',
    'Guard': 'Cotyledon/leaf (stomata)',
}

# ============================================================
# 1. CCC strength matrices for Ca2+ and K+
# ============================================================
ca2_matrix = ca2_comm.groupby(['Source', 'Target'])['Prob'].sum().reset_index()
ca2_pivot = ca2_matrix.pivot(index='Source', columns='Target', values='Prob').fillna(0)

k_matrix = k_comm.groupby(['Source', 'Target'])['Prob'].sum().reset_index()
k_pivot = k_matrix.pivot(index='Source', columns='Target', values='Prob').fillna(0)

ct_order = ['Epidermal', 'Mesophyll', 'Stele', 'Meristematic', 'Guard']
ca2_pivot = ca2_pivot.reindex(index=ct_order, columns=ct_order).fillna(0)
k_pivot = k_pivot.reindex(index=ct_order, columns=ct_order).fillna(0)

print('Ca2+ CCC strength matrix:')
print(ca2_pivot.round(2).to_string())
print()
print('K+ CCC strength matrix:')
print(k_pivot.round(2).to_string())

# ============================================================
# 2. Detect circuit: strongest Ca2+/K+ signaling pairs
# ============================================================
circuit_pairs = []
for src in ct_order:
    for tgt in ct_order:
        ca2_val = ca2_pivot.loc[src, tgt] if src in ca2_pivot.index and tgt in ca2_pivot.columns else 0
        k_val = k_pivot.loc[src, tgt] if src in k_pivot.index and tgt in k_pivot.columns else 0
        total = ca2_val + k_val
        if total > 0:
            circuit_pairs.append({
                'source': src, 'target': tgt,
                'source_organ': ct_to_organ.get(src, ''),
                'target_organ': ct_to_organ.get(tgt, ''),
                'ca2_strength': ca2_val, 'k_strength': k_val,
                'total_strength': total,
                'ca2_dominant': ca2_val > k_val
            })

circuit_df = pd.DataFrame(circuit_pairs).sort_values('total_strength', ascending=False)
circuit_df.to_csv(RESULTS + '/tables/ca2_k_ccc_circuit.csv', index=False)
print('\nSaved: ca2_k_ccc_circuit.csv')
print('\nTop 10 circuit pairs:')
print(circuit_df.head(10)[['source', 'target', 'ca2_strength', 'k_strength', 'total_strength']].to_string())

# ============================================================
# 3. Key Ca2+/K+ genes per cell type (for circuit annotation)
# ============================================================
expr = pd.read_csv(ATLAS + '/ca2_k_expr_per_celltype.csv')
expr_ct = expr[expr['celltype'] != 'ALL']

# Key circuit genes
circuit_genes = ['CBL9', 'CIPK23', 'AKT1', 'CML24', 'KC1', 'CDPK3', 'CBL1', 'CBL2', 'CIPK1', 'CIPK9']
gene_expr_per_ct = expr_ct[expr_ct['alias'].isin(circuit_genes)].pivot(
    index='alias', columns='celltype', values='mean_expr').reindex(circuit_genes)[ct_order]

print('\nKey circuit gene expression per cell type:')
print(gene_expr_per_ct.round(4).to_string())

# ============================================================
# 4. CCC Circuit Diagram
# ============================================================
print('\n=== Creating CCC circuit diagram ===')

fig, ax = plt.subplots(figsize=(12, 9))
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_aspect('equal')
ax.axis('off')

# Cell type positions (circular layout)
n_ct = len(ct_order)
angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n_ct, endpoint=False)
ct_positions = {}
ct_colors = {
    'Epidermal': '#0072B2', 'Mesophyll': '#009E73', 'Stele': '#D55E00',
    'Meristematic': '#CC79A7', 'Guard': '#E69F00'
}
radius = 3.5
for i, ct_name in enumerate(ct_order):
    ct_positions[ct_name] = (radius * np.cos(angles[i]), radius * np.sin(angles[i]))

# Draw edges (Ca2+ and K+)
max_strength = circuit_df['total_strength'].max()
for _, row in circuit_df.iterrows():
    src = row['source']
    tgt = row['target']
    if src not in ct_positions or tgt not in ct_positions:
        continue
    x1, y1 = ct_positions[src]
    x2, y2 = ct_positions[tgt]
    # Skip self-loops for clarity (draw as small arc)
    if src == tgt:
        continue

    ca2_val = row['ca2_strength']
    k_val = row['k_strength']

    # Draw Ca2+ edge (blue)
    if ca2_val > 0:
        width = 0.5 + 4 * (ca2_val / max_strength)
        alpha = 0.2 + 0.6 * (ca2_val / max_strength)
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#0072B2',
                                   lw=width, alpha=alpha, connectionstyle='arc3,rad=0.1'))

    # Draw K+ edge (green, offset)
    if k_val > 0:
        width = 0.3 + 3 * (k_val / max_strength)
        alpha = 0.2 + 0.6 * (k_val / max_strength)
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#009E73',
                                   lw=width, alpha=alpha, connectionstyle='arc3,rad=-0.15'))

# Draw cell type nodes
for ct_name, (x, y) in ct_positions.items():
    circle = plt.Circle((x, y), 0.6, color=ct_colors[ct_name], alpha=0.8, ec='black', linewidth=1.5)
    ax.add_patch(circle)
    ax.text(x, y, ct_name, ha='center', va='center', fontsize=8, fontweight='bold', color='white')
    # Organ proxy below
    ax.text(x, y - 0.9, ct_to_organ[ct_name], ha='center', va='top', fontsize=6,
            style='italic', color='0.4', wrap=True)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#0072B2', alpha=0.7, label='Ca2+ signaling'),
    mpatches.Patch(facecolor='#009E73', alpha=0.7, label='K+ signaling'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

ax.set_title('Ca2+/K+ Cell-Cell Communication Circuit\n(PlantCellChat, edge width = strength)',
             fontweight='bold', fontsize=13, pad=20)

fig.tight_layout()
fig.savefig(f'{FIGDIR}/ca2_k_ccc_circuit.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIGDIR}/ca2_k_ccc_circuit.svg', bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: ca2_k_ccc_circuit.png/.svg')

# ============================================================
# 5. Summary: the detected circuit
# ============================================================
print('\n=== DETECTED Ca2+/K+ SIGNALING CIRCUIT ===')
print('Top 5 Ca2+ circuits:')
ca2_top = circuit_df.nlargest(5, 'ca2_strength')
for _, r in ca2_top.iterrows():
    print(f"  {r['source']} -> {r['target']}: Ca2+={r['ca2_strength']:.2f}, K+={r['k_strength']:.2f}")

print('\nTop 5 K+ circuits:')
k_top = circuit_df.nlargest(5, 'k_strength')
for _, r in k_top.iterrows():
    print(f"  {r['source']} -> {r['target']}: Ca2+={r['ca2_strength']:.2f}, K+={r['k_strength']:.2f}")

print('\nCrosstalk (Ca2+ and K+ both active):')
both = circuit_df[(circuit_df['ca2_strength'] > 0) & (circuit_df['k_strength'] > 0)].sort_values('total_strength', ascending=False)
for _, r in both.head(5).iterrows():
    print(f"  {r['source']} -> {r['target']}: Ca2+={r['ca2_strength']:.2f}, K+={r['k_strength']:.2f}")

print('\nDone.')
