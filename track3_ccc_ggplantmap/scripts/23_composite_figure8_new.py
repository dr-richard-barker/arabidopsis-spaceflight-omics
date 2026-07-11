"""
Composite Figure 8 (new):
  Panel A (left):
    a1 (top-left):  ggpathway cascade network (CBL9-CIPK23-AKT1)
    a2 (bottom-left): ggPlantmap leaf cross-section (Ca2+/K+ circuit)
  Panel B (right): CBL9-CIPK23-AKT1 molecular cascade Sankey (matplotlib, no label clipping)

Layout:
  - Panel A: left 3000 px wide, full height 4640
  - Panel B: matplotlib Sankey scaled to fill right side
  - Panel labels: bold 'a' top-left, 'b' top-right
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

import subprocess, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# ============================================================
# 0. Dimensions
# ============================================================
CANVAS_H  = 4640
PANEL_A_W = 3000
MARGIN    = 40
LABEL_FONT_SIZE = 80

# ============================================================
# 1. Generate Panel B (matplotlib Sankey — no label clipping)
# ============================================================
result = subprocess.run(
    ['python3', ATLAS + '/make_sankey_panel_b.py'],
    capture_output=True, text=True
)
print(result.stdout[-500:] if result.stdout else '')
if result.returncode != 0:
    print('STDERR:', result.stderr[-500:])
    raise RuntimeError('Sankey generation failed')

sankey_path = WORK + '/figure8b_matplotlib.png'
print(f'Panel B: {sankey_path} ({os.path.getsize(sankey_path):,} bytes)')

# ============================================================
# 2. Load all panels
# ============================================================
a1_img = Image.open(WORK + '/figure8a1_cascade.png').convert('RGB')
a2_img = Image.open(WORK + '/figure8a2_leaf.png').convert('RGB')
b_img  = Image.open(sankey_path).convert('RGB')

print(f'a1 original: {a1_img.size}')
print(f'a2 original: {a2_img.size}')
print(f'b  original: {b_img.size}')

# ============================================================
# 3. Scale panel A sub-panels to fit 3000 wide x 4640 tall
# ============================================================
a1_w = PANEL_A_W
a1_h = int(a1_img.height * PANEL_A_W / a1_img.width)
a1_scaled = a1_img.resize((a1_w, a1_h), Image.LANCZOS)

a2_w = PANEL_A_W
a2_h = int(a2_img.height * PANEL_A_W / a2_img.width)
a2_scaled = a2_img.resize((a2_w, a2_h), Image.LANCZOS)

panel_a_content_h = a1_h + a2_h
print(f'a1 scaled: {a1_scaled.size}')
print(f'a2 scaled: {a2_scaled.size}')
print(f'Panel A content height: {panel_a_content_h} (canvas: {CANVAS_H})')

# Scale down if content exceeds canvas
if panel_a_content_h > CANVAS_H - 2 * MARGIN:
    sf = (CANVAS_H - 2 * MARGIN) / panel_a_content_h
    a1_scaled = a1_scaled.resize((PANEL_A_W, int(a1_h * sf)), Image.LANCZOS)
    a2_scaled = a2_scaled.resize((PANEL_A_W, int(a2_h * sf)), Image.LANCZOS)
    print(f'Scaled down: a1={a1_scaled.size}, a2={a2_scaled.size}')

# ============================================================
# 4. Scale panel B to canvas height (maintain aspect ratio)
# ============================================================
b_scale_h = CANVAS_H
b_scale_w = int(b_img.width * CANVAS_H / b_img.height)
b_scaled = b_img.resize((b_scale_w, b_scale_h), Image.LANCZOS)
print(f'Panel B scaled to: {b_scaled.size}')

CANVAS_W = PANEL_A_W + b_scale_w

# ============================================================
# 5. Composite onto canvas
# ============================================================
canvas = Image.new('RGB', (CANVAS_W, CANVAS_H), 'white')

# Panel A: a1 top, a2 bottom
canvas.paste(a1_scaled, (0, MARGIN))
a2_y = MARGIN + a1_scaled.height + MARGIN // 2
canvas.paste(a2_scaled, (0, a2_y))

# Panel B: right side
canvas.paste(b_scaled, (PANEL_A_W, 0))

# ============================================================
# 6. Panel labels and divider
# ============================================================
draw = ImageDraw.Draw(canvas)
try:
    font_label = ImageFont.truetype(
        '/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf', LABEL_FONT_SIZE)
except Exception:
    font_label = ImageFont.load_default()

draw.text((15, 10), 'a', fill='black', font=font_label)
draw.text((PANEL_A_W + 15, 10), 'b', fill='black', font=font_label)
draw.line([(PANEL_A_W, 0), (PANEL_A_W, CANVAS_H)], fill='#CCCCCC', width=3)

print(f'Final canvas size: {canvas.size}')

# ============================================================
# 7. Save
# ============================================================
out_png = WORK + '/Figure8_sankey_ca2_k_cascade_new.png'
canvas.save(out_png, dpi=(300, 300))
print(f'Saved: {out_png} ({os.path.getsize(out_png):,} bytes)')
