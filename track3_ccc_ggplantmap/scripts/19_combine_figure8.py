"""
Combine updated Figure 8:
  Panel A (left): ggKEGG CBL9-CIPK23-AKT1 cascade on ath04075 (from R script)
  Panel B (right): CBL9-CIPK23-AKT1 molecular cascade Sankey (regenerated standalone)

Approach: regenerate panel B alone via plotly, then composite with panel A (PIL).
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
import pandas as pd
import numpy as np
import os
from PIL import Image, ImageDraw, ImageFont

# ============================================================
# Regenerate Panel B (molecular cascade Sankey) as standalone
# ============================================================
expr = pd.read_csv(ATLAS + '/ca2_k_expr_per_celltype.csv')
coexpr_guard = pd.read_csv(RESULTS + '/tables/ca2_k_coexpression_Guard.csv', index_col=0)
coexpr_guard.index = [x.replace('_Guard', '') for x in coexpr_guard.index]


def get_guard_expr(gene_alias):
    sub = expr[(expr['alias'] == gene_alias) & (expr['celltype'] == 'Guard')]
    if len(sub) > 0:
        return sub.iloc[0]['mean_expr']
    return 0.01


def scale_expr(val):
    return max(val * 100, 1.0)


def get_coexpr(g1, g2):
    if g1 in coexpr_guard.index and g2 in coexpr_guard.columns:
        return coexpr_guard.loc[g1, g2]
    return 0


cascade_nodes = [
    'Ca2+', 'CBL1', 'CBL9', 'CaM3', 'CML24', 'CIPK23',
    'CDPK4', 'CDPK3', 'AKT1', 'ANNAT1', 'KC1', 'K+ uptake',
]
node_colors_b = [
    '#E9ED4C', '#0279EE', '#0279EE', '#FD9BED', '#FD9BED', '#75A025',
    '#75A025', '#75A025', '#FF9400', '#FF9400', '#FF9400', '#ECE9E2',
]

cascade_links = [
    {'source': 0, 'target': 1, 'value': scale_expr(get_guard_expr('CBL1')),
     'color': 'rgba(233, 237, 76, 0.5)'},
    {'source': 0, 'target': 2, 'value': scale_expr(get_guard_expr('CBL9')),
     'color': 'rgba(233, 237, 76, 0.5)'},
    {'source': 0, 'target': 3, 'value': scale_expr(get_guard_expr('CaM3')),
     'color': 'rgba(233, 237, 76, 0.5)'},
    {'source': 0, 'target': 4, 'value': scale_expr(get_guard_expr('CML24')),
     'color': 'rgba(233, 237, 76, 0.5)'},
    {'source': 1, 'target': 5, 'value': scale_expr(get_guard_expr('CIPK23')),
     'color': 'rgba(2, 121, 238, 0.4)'},
    {'source': 2, 'target': 5, 'value': scale_expr(get_guard_expr('CIPK23')),
     'color': 'rgba(2, 121, 238, 0.4)'},
    {'source': 3, 'target': 6, 'value': scale_expr(get_guard_expr('CDPK4')),
     'color': 'rgba(253, 155, 237, 0.4)'},
    {'source': 3, 'target': 7, 'value': scale_expr(get_guard_expr('CDPK3')),
     'color': 'rgba(253, 155, 237, 0.4)'},
    {'source': 5, 'target': 8, 'value': scale_expr(get_guard_expr('AKT1')),
     'color': 'rgba(255, 148, 0, 0.6)'},
    {'source': 6, 'target': 9, 'value': scale_expr(get_guard_expr('ANNAT1')),
     'color': 'rgba(117, 160, 37, 0.4)'},
    {'source': 4, 'target': 10, 'value': scale_expr(get_guard_expr('KC1')),
     'color': 'rgba(253, 155, 237, 0.4)'},
    {'source': 8, 'target': 11, 'value': scale_expr(get_guard_expr('AKT1')),
     'color': 'rgba(255, 148, 0, 0.6)'},
    {'source': 10, 'target': 11, 'value': scale_expr(get_guard_expr('KC1')),
     'color': 'rgba(255, 148, 0, 0.4)'},
]

# Edge labels showing co-expression r for key edges
link_labels = [
    '', '', '', '',
    f"r={get_coexpr('CBL1','CIPK23'):.2f}",
    f"r={get_coexpr('CBL9','CIPK23'):.2f}",
    f"r={get_coexpr('CaM3','CDPK4'):.2f}",
    f"r={get_coexpr('CaM3','CDPK3'):.2f}",
    f"r={get_coexpr('CIPK23','AKT1'):.2f}",
    f"r={get_coexpr('CDPK4','ANNAT1'):.2f}",
    f"r={get_coexpr('CML24','KC1'):.2f}",
    '', '',
]

sankey_b = go.Sankey(
    node=dict(
        pad=30, thickness=22,
        line=dict(color='black', width=0.5),
        label=cascade_nodes, color=node_colors_b,
    ),
    link=dict(
        source=[l['source'] for l in cascade_links],
        target=[l['target'] for l in cascade_links],
        value=[l['value'] for l in cascade_links],
        color=[l['color'] for l in cascade_links],
        label=link_labels,
        hovertemplate='%{label}<extra></extra>',
    ),
    arrangement='snap',
)

fig_b = go.Figure(data=[sankey_b])
fig_b.update_layout(
    font=dict(family='Arial, sans-serif', size=16),
    width=1000, height=1500,
    margin=dict(l=20, r=20, t=20, b=20),
)

panel_b_file = WORK + '/figure8b_sankey.png'
fig_b.write_image(panel_b_file, width=1000, height=1500, scale=3)
print(f"Panel B saved: {os.path.getsize(panel_b_file)/1024:.0f} KB")

# ============================================================
# Composite panel A + panel B
# ============================================================
panel_a = Image.open(WORK + '/figure8a_ggkegg.png').convert('RGB')
panel_b = Image.open(panel_b_file).convert('RGB')

# Target height: match to a common value
target_h = 4500
# Scale panel A (already 3000x4500)
a_ratio = target_h / panel_a.height
panel_a_r = panel_a.resize((int(panel_a.width * a_ratio), target_h), Image.LANCZOS)
# Scale panel B (3000x4500 at scale=3 from 1000x1500)
b_ratio = target_h / panel_b.height
panel_b_r = panel_b.resize((int(panel_b.width * b_ratio), target_h), Image.LANCZOS)

# Gap between panels
gap = 80
label_h = 140  # space for panel labels at top
total_w = panel_a_r.width + gap + panel_b_r.width
total_h = target_h + label_h

composite = Image.new('RGB', (total_w, total_h), 'white')
composite.paste(panel_a_r, (0, label_h))
composite.paste(panel_b_r, (panel_a_r.width + gap, label_h))

# Add panel labels
draw = ImageDraw.Draw(composite)
try:
    font_lbl = ImageFont.truetype('/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf', 90)
    font_sub = ImageFont.truetype('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf', 48)
except Exception:
    font_lbl = ImageFont.load_default()
    font_sub = ImageFont.load_default()

# Panel a label
draw.text((30, 20), 'a', fill='black', font=font_lbl)
draw.text((130, 55), 'KEGG pathway overlay (ath04075)', fill='#444444', font=font_sub)
# Panel b label
bx = panel_a_r.width + gap
draw.text((bx + 30, 20), 'b', fill='black', font=font_lbl)
draw.text((bx + 130, 55), 'Molecular cascade (guard cell expression)', fill='#444444', font=font_sub)

# Save composite
out_png = RESULTS + '/figures/Figure8_sankey_ca2_k_cascade.png'
composite.save(out_png, dpi=(300, 300))

# Verify
arr = np.array(composite)
nw = (arr < 250).any(axis=2).mean() * 100
print(f"Composite Figure 8: {composite.size}, {os.path.getsize(out_png)/1024:.0f} KB, non-white={nw:.1f}%")
print("Done.")
