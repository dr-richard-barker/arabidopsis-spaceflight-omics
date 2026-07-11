"""
Step 1: Two-panel Sankey diagram (Figure 8).
Panel A: Cell-type Ca2+/K+ signal flow (inter-cell-type only, self-loops excluded)
Panel B: CBL9-CIPK23-AKT1 molecular cascade with Guard cell expression values
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

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

# ============================================================
# Load data
# ============================================================
ccc = pd.read_csv(RESULTS + '/tables/ccc_per_signal_strength.csv')
expr = pd.read_csv(ATLAS + '/ca2_k_expr_per_celltype.csv')
coexpr_guard = pd.read_csv(RESULTS + '/tables/ca2_k_coexpression_Guard.csv', index_col=0)
coexpr_guard.index = [x.replace('_Guard', '') for x in coexpr_guard.index]

# ============================================================
# Panel A: Cell-type signal flow Sankey
# ============================================================
# Filter Ca2+ and K+ inter-cell-type flows (exclude self-loops)
ca2_flows = ccc[(ccc['Signal'] == 'Ca2+') & (ccc['Source'] != ccc['Target'])].copy()
k_flows = ccc[(ccc['Signal'] == 'K+') & (ccc['Source'] != ccc['Target'])].copy()

cell_types = ['Guard', 'Meristematic', 'Epidermal', 'Mesophyll', 'Stele']
# Phylo palette colors
ct_colors = {
    'Guard': '#FF9400',     # orange (stomatal)
    'Meristematic': '#75A025',  # green
    'Epidermal': '#0279EE',     # blue
    'Mesophyll': '#FD9BED',     # pink
    'Stele': '#E9ED4C',         # yellow
}

# Sankey nodes: sources (0-4) + targets (5-9)
node_labels_a = [f'{ct} (src)' for ct in cell_types] + [f'{ct} (tgt)' for ct in cell_types]
node_colors_a = [ct_colors[ct] for ct in cell_types] + [ct_colors[ct] for ct in cell_types]

# Build links: Ca2+ flows (blue) then K+ flows (green)
links_a = []
for _, r in ca2_flows.iterrows():
    src_idx = cell_types.index(r['Source'])
    tgt_idx = cell_types.index(r['Target']) + 5
    links_a.append({'source': src_idx, 'target': tgt_idx, 'value': r['Prob'],
                    'color': 'rgba(2, 121, 238, 0.4)',  # blue transparent
                    'label': f"Ca2+: {r['Source']}→{r['Target']} ({r['Prob']:.1f})"})

for _, r in k_flows.iterrows():
    src_idx = cell_types.index(r['Source'])
    tgt_idx = cell_types.index(r['Target']) + 5
    links_a.append({'source': src_idx, 'target': tgt_idx, 'value': r['Prob'],
                    'color': 'rgba(117, 160, 37, 0.4)',  # green transparent
                    'label': f"K+: {r['Source']}→{r['Target']} ({r['Prob']:.2f})"})

# Scale K+ values for visibility (they're 10x smaller)
# Use a minimum width and log-ish scaling
for link in links_a:
    if 'K+' in link['label']:
        link['value'] = max(link['value'], 0.3)  # minimum visibility

sankey_a = go.Sankey(
    node=dict(
        pad=25,
        thickness=20,
        line=dict(color='black', width=0.5),
        label=node_labels_a,
        color=node_colors_a
    ),
    link=dict(
        source=[l['source'] for l in links_a],
        target=[l['target'] for l in links_a],
        value=[l['value'] for l in links_a],
        color=[l['color'] for l in links_a],
        label=[l['label'] for l in links_a],
        hovertemplate='%{label}<extra></extra>'
    ),
    arrangement='snap',
    domain=dict(x=[0, 0.48], y=[0, 1])
)

# ============================================================
# Panel B: Molecular cascade Sankey
# ============================================================
# Get Guard cell expression for cascade genes
def get_guard_expr(gene_alias):
    sub = expr[(expr['alias'] == gene_alias) & (expr['celltype'] == 'Guard')]
    if len(sub) > 0:
        return sub.iloc[0]['mean_expr']
    return 0.01  # minimum for visibility

# Cascade structure:
# Ca2+ → CBL1, CBL9 (sensors)
# Ca2+ → CaM3, CML24 (calmodulin sensors)
# CBL1/CBL9 → CIPK23 (kinase)
# CaM3 → CDPK4, CDPK3 (kinases) [via Ca2+/CaM activation]
# CIPK23 → AKT1 (K+ channel)
# CDPK4 → ANNAT1 (annexin/Ca2+ buffering)
# AKT1 → K+ uptake
# CML24 → KC1 (regulation)

# Node definitions
cascade_nodes = [
    'Ca2+',           # 0
    'CBL1',           # 1
    'CBL9',           # 2
    'CaM3',           # 3
    'CML24',          # 4
    'CIPK23',         # 5
    'CDPK4',          # 6
    'CDPK3',          # 7
    'AKT1',           # 8
    'ANNAT1',         # 9
    'KC1',            # 10
    'K+ uptake',      # 11
]

node_colors_b = [
    '#E9ED4C',  # Ca2+ - yellow
    '#0279EE',  # CBL1 - blue
    '#0279EE',  # CBL9 - blue
    '#FD9BED',  # CaM3 - pink
    '#FD9BED',  # CML24 - pink
    '#75A025',  # CIPK23 - green
    '#75A025',  # CDPK4 - green
    '#75A025',  # CDPK3 - green
    '#FF9400',  # AKT1 - orange
    '#FF9400',  # ANNAT1 - orange
    '#FF9400',  # KC1 - orange
    '#ECE9E2',  # K+ uptake - light
]

# Build links with expression values as flow widths
# Scale expression for visibility (multiply by 100, min 1)
def scale_expr(val):
    return max(val * 100, 1.0)

# Get co-expression r values for labels
def get_coexpr(g1, g2):
    if g1 in coexpr_guard.index and g2 in coexpr_guard.columns:
        return coexpr_guard.loc[g1, g2]
    return 0

cascade_links = [
    # Ca2+ → sensors
    {'source': 0, 'target': 1, 'value': scale_expr(get_guard_expr('CBL1')),
     'color': 'rgba(233, 237, 76, 0.5)', 'label': f"Ca2+→CBL1 (expr={get_guard_expr('CBL1'):.3f})"},
    {'source': 0, 'target': 2, 'value': scale_expr(get_guard_expr('CBL9')),
     'color': 'rgba(233, 237, 76, 0.5)', 'label': f"Ca2+→CBL9 (expr={get_guard_expr('CBL9'):.3f})"},
    {'source': 0, 'target': 3, 'value': scale_expr(get_guard_expr('CaM3')),
     'color': 'rgba(233, 237, 76, 0.5)', 'label': f"Ca2+→CaM3 (expr={get_guard_expr('CaM3'):.3f})"},
    {'source': 0, 'target': 4, 'value': scale_expr(get_guard_expr('CML24')),
     'color': 'rgba(233, 237, 76, 0.5)', 'label': f"Ca2+→CML24 (expr={get_guard_expr('CML24'):.3f})"},

    # CBL1/CBL9 → CIPK23
    {'source': 1, 'target': 5, 'value': scale_expr(get_guard_expr('CIPK23')),
     'color': 'rgba(2, 121, 238, 0.4)', 'label': f"CBL1→CIPK23 (r={get_coexpr('CBL1','CIPK23'):.3f})"},
    {'source': 2, 'target': 5, 'value': scale_expr(get_guard_expr('CIPK23')),
     'color': 'rgba(2, 121, 238, 0.4)', 'label': f"CBL9→CIPK23 (r={get_coexpr('CBL9','CIPK23'):.3f})"},

    # CaM3 → CDPK4, CDPK3
    {'source': 3, 'target': 6, 'value': scale_expr(get_guard_expr('CDPK4')),
     'color': 'rgba(253, 155, 237, 0.4)', 'label': f"CaM3→CDPK4 (r={get_coexpr('CaM3','CDPK4'):.3f})"},
    {'source': 3, 'target': 7, 'value': scale_expr(get_guard_expr('CDPK3')),
     'color': 'rgba(253, 155, 237, 0.4)', 'label': f"CaM3→CDPK3 (r={get_coexpr('CaM3','CDPK3'):.3f})"},

    # CIPK23 → AKT1 (key crosstalk edge)
    {'source': 5, 'target': 8, 'value': scale_expr(get_guard_expr('AKT1')),
     'color': 'rgba(255, 148, 0, 0.6)', 'label': f"CIPK23→AKT1 (r={get_coexpr('CIPK23','AKT1'):.3f})"},

    # CDPK4 → ANNAT1
    {'source': 6, 'target': 9, 'value': scale_expr(get_guard_expr('ANNAT1')),
     'color': 'rgba(117, 160, 37, 0.4)', 'label': f"CDPK4→ANNAT1 (r={get_coexpr('CDPK4','ANNAT1'):.3f})"},

    # CML24 → KC1
    {'source': 4, 'target': 10, 'value': scale_expr(get_guard_expr('KC1')),
     'color': 'rgba(253, 155, 237, 0.4)', 'label': f"CML24→KC1 (r={get_coexpr('CML24','KC1'):.3f})"},

    # AKT1 → K+ uptake
    {'source': 8, 'target': 11, 'value': scale_expr(get_guard_expr('AKT1')),
     'color': 'rgba(255, 148, 0, 0.6)', 'label': "AKT1→K+ uptake"},
    # KC1 → K+ uptake (KC1 modulates AKT1)
    {'source': 10, 'target': 11, 'value': scale_expr(get_guard_expr('KC1')),
     'color': 'rgba(255, 148, 0, 0.4)', 'label': "KC1→K+ uptake"},
]

sankey_b = go.Sankey(
    node=dict(
        pad=25,
        thickness=20,
        line=dict(color='black', width=0.5),
        label=cascade_nodes,
        color=node_colors_b
    ),
    link=dict(
        source=[l['source'] for l in cascade_links],
        target=[l['target'] for l in cascade_links],
        value=[l['value'] for l in cascade_links],
        color=[l['color'] for l in cascade_links],
        label=[l['label'] for l in cascade_links],
        hovertemplate='%{label}<extra></extra>'
    ),
    arrangement='snap',
    domain=dict(x=[0.52, 1.0], y=[0, 1])
)

# ============================================================
# Combine into single figure
# ============================================================
fig = go.Figure(data=[sankey_a, sankey_b])

fig.update_layout(
    title=dict(
        text="<b>Ca2+/K+ signaling circuit: signal flow and molecular cascade</b>",
        x=0.5, font=dict(size=18, family='Arial, sans-serif')
    ),
    font=dict(family='Arial, sans-serif', size=12),
    width=1800,
    height=900,
    margin=dict(l=20, r=20, t=60, b=20),
    annotations=[
        dict(text="<b>a) Cell-type signal flow</b><br>(Ca2+ blue, K+ green; width = CCC strength)",
             x=0.24, y=-0.02, xref='paper', yref='paper',
             showarrow=False, font=dict(size=13)),
        dict(text="<b>b) CBL9-CIPK23-AKT1 cascade</b><br>(guard cell expression; labels show co-expression r)",
             x=0.76, y=-0.02, xref='paper', yref='paper',
             showarrow=False, font=dict(size=13)),
    ]
)

# Export
out_png = RESULTS + '/figures/Figure8_sankey_ca2_k_cascade.png'
out_svg = RESULTS + '/figures/Figure8_sankey_ca2_k_cascade.svg'
fig.write_image(out_png, width=1800, height=900, scale=3)
fig.write_image(out_svg, width=1800, height=900)

import os
from PIL import Image
img = Image.open(out_png)
arr = np.array(img)
nw = (arr < 250).any(axis=2).sum()
print(f"Figure 8 Sankey: {img.size}, {os.path.getsize(out_png)/1024:.0f} KB, non-white={100*nw/(arr.shape[0]*arr.shape[1]):.1f}%")
print(f"SVG: {os.path.getsize(out_svg)/1024:.0f} KB")
print("Done.")
