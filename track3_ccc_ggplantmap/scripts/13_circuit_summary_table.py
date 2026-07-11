"""
Step 6: Ca2+/K+ signaling circuit summary table.
Synthesizes CCC strength, co-expression, and expression data into a single
circuit-level summary describing the cell-to-cell signaling circuit connected
to the Ca2+/K+ pathway diagram.
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

import pandas as pd
import numpy as np

# --- Load data ---
ccc = pd.read_csv(RESULTS + '/tables/ca2_k_ccc_circuit.csv')
coexpr = pd.read_csv(RESULTS + '/tables/ca2_k_coexpression_matrix.csv', index_col=0)
expr = pd.read_csv(ATLAS + '/ca2_k_expr_per_celltype.csv')

# Per-cell-type co-expression matrices
coexpr_by_ct = {}
for ct in ['Epidermal', 'Mesophyll', 'Stele', 'Meristematic', 'Guard']:
    coexpr_by_ct[ct] = pd.read_csv(RESULTS + f'/tables/ca2_k_coexpression_{ct}.csv', index_col=0)

# Expression pivot: alias -> celltype -> (mean, pct)
expr_pivot = expr.pivot_table(index='alias', columns='celltype', values='mean_expr')
pct_pivot = expr.pivot_table(index='alias', columns='celltype', values='pct_expr')

# --- Key circuit genes per cell type (top 3 by expression) ---
circuit_genes = ['CBL1','CBL9','CBL2','CIPK23','CIPK1','CIPK9','CML38','CML42',
                 'CML24','CaM7','CaM3','CDPK3','CDPK4','CAX3','NRT2.1','SOS1',
                 'ANNAT1','AKT1','AKT2','GORK','KAT1','KAT2','KC1']

def top_genes_for_ct(ct, n=3):
    sub = expr[(expr['celltype']==ct) & (expr['alias'].isin(circuit_genes))]
    sub = sub.sort_values('mean_expr', ascending=False)
    return ', '.join(f"{r['alias']}({r['mean_expr']:.2f})" for _, r in sub.head(n).iterrows())

# --- Build circuit summary rows ---
# For each CCC pair (source -> target), summarize the circuit
rows = []
for _, r in ccc.iterrows():
    src = r['source']
    tgt = r['target']
    ca2_s = r['ca2_strength']
    k_s = r['k_strength']
    total = r['total_strength']

    # Key ligand genes in source (highest expression)
    src_genes = top_genes_for_ct(src, 3)
    tgt_genes = top_genes_for_ct(tgt, 3)

    # Strongest co-expression pair in source cell type
    src_coexpr = coexpr_by_ct.get(src)
    if src_coexpr is not None:
        # Get upper triangle max
        mask = np.triu(np.ones_like(src_coexpr, dtype=bool), k=1)
        vals = src_coexpr.where(mask).stack()
        if len(vals) > 0:
            top_pair = vals.idxmax()
            top_r = vals.max()
            coexpr_str = f"{top_pair[0]}~{top_pair[1]} (r={top_r:.3f})"
        else:
            coexpr_str = "N/A"
    else:
        coexpr_str = "N/A"

    # Dominant signal
    dominant = "Ca2+" if ca2_s > k_s else "K+"
    ratio = ca2_s / k_s if k_s > 0 else float('inf')

    rows.append({
        'source_cell_type': src,
        'target_cell_type': tgt,
        'source_organ_proxy': r['source_organ'],
        'target_organ_proxy': r['target_organ'],
        'dominant_signal': dominant,
        'ca2_strength': round(ca2_s, 2),
        'k_strength': round(k_s, 2),
        'total_strength': round(total, 2),
        'ca2_k_ratio': round(ratio, 1),
        'top_source_genes': src_genes,
        'top_target_genes': tgt_genes,
        'strongest_coexpr_in_source': coexpr_str,
    })

summary = pd.DataFrame(rows)
summary = summary.sort_values('total_strength', ascending=False)
summary.to_csv(RESULTS + '/tables/ca2_k_signaling_circuit.csv', index=False)
print(f"Circuit summary: {summary.shape[0]} pairs")
print(summary[['source_cell_type','target_cell_type','dominant_signal',
               'ca2_strength','k_strength','total_strength','strongest_coexpr_in_source']].to_string(index=False))

# --- Also create a gene-level circuit node table ---
node_rows = []
for g in circuit_genes:
    row = {'gene': g}
    for ct in ['Epidermal','Mesophyll','Stele','Meristematic','Guard']:
        row[f'{ct}_mean'] = round(expr_pivot.loc[g, ct], 3) if g in expr_pivot.index else 0
        row[f'{ct}_pct'] = round(pct_pivot.loc[g, ct], 1) if g in pct_pivot.index else 0
    # Max co-expression partner globally
    if g in coexpr.index:
        others = coexpr.loc[g].drop(g)
        top_partner = others.idxmax()
        row['top_coexpr_partner'] = top_partner
        row['top_coexpr_r'] = round(others.max(), 3)
    else:
        row['top_coexpr_partner'] = ''
        row['top_coexpr_r'] = 0
    node_rows.append(row)

nodes = pd.DataFrame(node_rows)
nodes.to_csv(RESULTS + '/tables/ca2_k_circuit_nodes.csv', index=False)
print(f"\nCircuit nodes: {nodes.shape[0]} genes")
print(nodes[['gene','Guard_mean','Guard_pct','top_coexpr_partner','top_coexpr_r']].to_string(index=False))
