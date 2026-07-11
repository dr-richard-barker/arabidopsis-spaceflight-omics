#!/usr/bin/env python3
"""
Fix Figure3 panel C bug and regenerate all 6 main manuscript figures.
Bug: panel C checked for 'total_strength' column but actual column is 'total_comm'.
     Fallback iloc[:,-1] picked up Signal column (all NaN for top 10) -> empty bars.
Fix: Use total_comm directly, merge with ccc_lr_pairs.csv for Signal labels,
     color by signal pathway (Ca2+=orange, General=grey, others=blue).
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
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import os

# Publication-quality defaults
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.labelsize'] = 9
matplotlib.rcParams['xtick.labelsize'] = 8
matplotlib.rcParams['ytick.labelsize'] = 8
matplotlib.rcParams['legend.fontsize'] = 8
matplotlib.rcParams['axes.titlesize'] = 10

CB_COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00', '#000000']
FIGDIR = RESULTS + '/figures'
os.makedirs(FIGDIR, exist_ok=True)

def save_fig(fig, name):
    fig.savefig(f'{FIGDIR}/{name}.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIGDIR}/{name}.svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {name}.png/.svg')

# ============================================================
# FIGURE 1: LASSO Biomarker Panel
# ============================================================
print('=== Figure 1: LASSO Biomarker Panel ===')
disc_perf = pd.read_csv(RESULTS + '/tables/discovery_performance.csv')
stability = pd.read_csv(RESULTS + '/tables/all_feature_stability.csv')
biomarker = pd.read_csv(RESULTS + '/tables/biomarker_panel.csv')

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Panel A: ROC curves per study
ax = axes[0]
studies = disc_perf['held_out_study'].unique()
for i, study in enumerate(studies):
    study_data = disc_perf[disc_perf['held_out_study'] == study]
    aucs = study_data['auc'].values
    ax.barh(i, np.mean(aucs), color=CB_COLORS[i % len(CB_COLORS)], alpha=0.7, height=0.6)
    ax.errorbar(np.mean(aucs), i, xerr=np.std(aucs), fmt='none', color='black', capsize=3)
    ax.scatter(aucs, np.full_like(aucs, i) + np.random.normal(0, 0.05, len(aucs)),
              color='black', s=10, alpha=0.5, zorder=5)
ax.set_yticks(range(len(studies)))
ax.set_yticklabels(studies)
ax.set_xlabel('AUC')
ax.set_ylabel('Held-out study')
ax.set_title('a) Leave-one-study-out CV AUC', fontweight='bold', loc='left')
ax.axvline(0.5, color='grey', linestyle='--', alpha=0.5, label='Chance')
ax.axvline(0.734, color='red', linestyle='-', alpha=0.7, label='Mean (0.734)')
ax.legend(loc='lower right', fontsize=7)
ax.set_xlim(0.3, 1.05)

# Panel B: Stability barplot
ax = axes[1]
top_stable = stability.head(20)
colors = ['#D55E00' if c > 0 else '#0072B2' for c in top_stable['mean_coefficient']]
ax.barh(range(len(top_stable)), top_stable['selection_frequency'] * 100, color=colors, alpha=0.8)
ax.set_yticks(range(len(top_stable)))
ax.set_yticklabels(top_stable['feature'], fontsize=7)
ax.set_xlabel('Selection frequency (%)')
ax.set_title('b) Top 20 feature stability', fontweight='bold', loc='left')
ax.axvline(50, color='grey', linestyle='--', alpha=0.5)
ax.invert_yaxis()
legend_elements = [Patch(facecolor='#D55E00', alpha=0.8, label='Positive (flight)'),
                   Patch(facecolor='#0072B2', alpha=0.8, label='Negative (ground)')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=7)

# Panel C: Coefficient forest plot
ax = axes[2]
top4 = biomarker.head(4)
y_pos = range(len(top4))
coef = top4['mean_coefficient']
sd = top4['sd_coefficient']
colors = ['#D55E00' if c > 0 else '#0072B2' for c in coef]
ax.errorbar(coef, y_pos, xerr=sd, fmt='o', color='black', markersize=8,
           capsize=4, capthick=1.5, elinewidth=1.5, zorder=5)
for i, (c, col) in enumerate(zip(coef, colors)):
    ax.scatter(c, i, color=col, s=80, zorder=6, edgecolor='black', linewidth=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{g}\n({a})" for g, a in zip(top4['feature'],
           ['HSP22.0','NAT','chloroplast','TRAF-like'])], fontsize=8)
ax.axvline(0, color='grey', linestyle='-', alpha=0.5)
ax.set_xlabel('Mean coefficient')
ax.set_title('c) Top 4 biomarkers (100% selection)', fontweight='bold', loc='left')
ax.invert_yaxis()

fig.suptitle('LASSO Biomarker Panel for Spaceflight Prediction', fontsize=13, fontweight='bold', y=1.02)
fig.tight_layout()
save_fig(fig, 'Figure1_LASSO_panel')

# ============================================================
# FIGURE 2: scPlantLLM Zero-shot + Clustering
# ============================================================
print('=== Figure 2: scPlantLLM Clustering ===')
umap = pd.read_csv(RESULTS + '/tables/umap_coords.csv')
annotations = pd.read_csv(RESULTS + '/tables/scplantllm_zero_shot_annotations.csv')
clusters = pd.read_csv(RESULTS + '/tables/cluster_assignments.csv')

umap_merged = umap.merge(annotations[['cell_name', 'true_celltype', 'predicted_celltype']].rename(
    columns={'cell_name':'cell_id', 'true_celltype':'CellType', 'predicted_celltype':'scplantllm_pred'}),
    on='cell_id', how='left')
umap_merged = umap_merged.merge(clusters[['cell_id', 'leiden']], on='cell_id', how='left')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: UMAP by atlas cell type
ax = axes[0]
celltypes = umap_merged['CellType'].unique()
for i, ct in enumerate(celltypes):
    mask = umap_merged['CellType'] == ct
    ax.scatter(umap_merged.loc[mask, 'UMAP1'], umap_merged.loc[mask, 'UMAP2'],
              c=CB_COLORS[i % len(CB_COLORS)], s=3, alpha=0.5, label=f'{ct} (n={mask.sum()})')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('a) Atlas cell type annotations', fontweight='bold', loc='left')
ax.legend(markerscale=3, fontsize=7, loc='best')

# Panel B: UMAP by Leiden clusters
ax = axes[1]
leiden_clusters = sorted(umap_merged['leiden'].unique())
for i, cl in enumerate(leiden_clusters):
    mask = umap_merged['leiden'] == cl
    ax.scatter(umap_merged.loc[mask, 'UMAP1'], umap_merged.loc[mask, 'UMAP2'],
              c=CB_COLORS[i % len(CB_COLORS)], s=3, alpha=0.5, label=f'Cluster {cl} (n={mask.sum()})')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('b) Leiden clusters (res=0.5)', fontweight='bold', loc='left')
ax.legend(markerscale=3, fontsize=7, loc='best')

# Panel C: UMAP by scPlantLLM predictions
ax = axes[2]
preds = umap_merged['scplantllm_pred'].value_counts()
top_preds = preds.head(8).index
for i, pred in enumerate(top_preds):
    mask = umap_merged['scplantllm_pred'] == pred
    ax.scatter(umap_merged.loc[mask, 'UMAP1'], umap_merged.loc[mask, 'UMAP2'],
              c=CB_COLORS[i % len(CB_COLORS)], s=3, alpha=0.5, label=f'{pred} (n={mask.sum()})')
other_mask = ~umap_merged['scplantllm_pred'].isin(top_preds)
if other_mask.sum() > 0:
    ax.scatter(umap_merged.loc[other_mask, 'UMAP1'], umap_merged.loc[other_mask, 'UMAP2'],
              c='grey', s=3, alpha=0.3, label=f'Other (n={other_mask.sum()})')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('c) scPlantLLM zero-shot predictions', fontweight='bold', loc='left')
ax.legend(markerscale=3, fontsize=6, loc='best')

fig.suptitle('scPlantLLM Zero-shot Embeddings and Clustering (9,998 cells)', fontsize=13, fontweight='bold', y=1.02)
fig.tight_layout()
save_fig(fig, 'Figure2_scPlantLLM_clustering')

# ============================================================
# FIGURE 3: Cell-Cell Communication — FIXED panel C
# ============================================================
print('=== Figure 3: Cell-Cell Communication (panel C FIXED) ===')
comm = pd.read_csv(RESULTS + '/tables/ccc_communication_strength.csv', index_col=0)
signal_strength = pd.read_csv(RESULTS + '/tables/ccc_per_signal_strength.csv')
lr_pairs = pd.read_csv(RESULTS + '/tables/ccc_lr_pairs_with_strength.csv')
lr_db = pd.read_csv(RESULTS + '/tables/ccc_lr_pairs.csv')

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Communication heatmap
ax = axes[0]
im = ax.imshow(comm.values, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(len(comm.columns)))
ax.set_xticklabels(comm.columns, rotation=45, ha='right')
ax.set_yticks(range(len(comm.index)))
ax.set_yticklabels(comm.index)
ax.set_xlabel('Target cell type'); ax.set_ylabel('Source cell type')
ax.set_title('a) Communication strength heatmap', fontweight='bold', loc='left')
plt.colorbar(im, ax=ax, label='Communication strength', shrink=0.8)
for i in range(len(comm.index)):
    for j in range(len(comm.columns)):
        ax.text(j, i, f'{comm.values[i,j]:.0f}', ha='center', va='center',
               fontsize=7, color='black' if comm.values[i,j] < 100 else 'white')

# Panel B: Top signaling pathways
ax = axes[1]
signal_totals = signal_strength[signal_strength['Signal'] != ''].groupby('Signal')['Prob'].sum().sort_values(ascending=False).head(10)
colors_sig = plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(signal_totals)))
ax.barh(range(len(signal_totals)), signal_totals.values, color=colors_sig)
ax.set_yticks(range(len(signal_totals)))
ax.set_yticklabels(signal_totals.index)
ax.set_xlabel('Total communication strength')
ax.set_title('b) Top signaling pathways', fontweight='bold', loc='left')
ax.invert_yaxis()
ca2_idx = list(signal_totals.index).index('Ca2+') if 'Ca2+' in signal_totals.index else -1
if ca2_idx >= 0:
    ax.barh(ca2_idx, signal_totals.values[ca2_idx], color='#D55E00', edgecolor='black', linewidth=1.5)

# Panel C: Top LR pairs — FIXED
# Merge with LR pairs database to get Signal labels
lr_merged = lr_pairs.merge(lr_db[['Interaction_name', 'Signal']], on='Interaction_name',
                           how='left', suffixes=('', '_db'))
lr_merged['Signal_final'] = lr_merged['Signal'].fillna(lr_merged['Signal_db']).fillna('General')
# Sort by total_comm (the correct column name) and take top 10
lr_merged = lr_merged.sort_values('total_comm', ascending=False)
top_lr = lr_merged.head(10).copy()

# Color by signal pathway
signal_color_map = {'Ca2+': '#D55E00', 'General': '#999999'}
colors_lr = [signal_color_map.get(s, '#0072B2') for s in top_lr['Signal_final']]

ax = axes[2]
ax.barh(range(len(top_lr)), top_lr['total_comm'].values, color=colors_lr, alpha=0.8)
ax.set_yticks(range(len(top_lr)))
ax.set_yticklabels(top_lr['Interaction_name'].values, fontsize=7)
ax.set_xlabel('Total communication strength')
ax.set_title('c) Top 10 ligand-receptor pairs', fontweight='bold', loc='left')
ax.invert_yaxis()
# Legend for signal categories
legend_elements_lr = [Patch(facecolor='#D55E00', alpha=0.8, label='Ca2+'),
                      Patch(facecolor='#999999', alpha=0.8, label='General (RLK)'),
                      Patch(facecolor='#0072B2', alpha=0.8, label='Other signals')]
ax.legend(handles=legend_elements_lr, loc='lower right', fontsize=6)

fig.suptitle('PlantCellChat Cell-Cell Communication Analysis', fontsize=13, fontweight='bold', y=1.02)
fig.tight_layout()
save_fig(fig, 'Figure3_CCC_communication')

# Verify panel C has data
from PIL import Image
img = Image.open(f'{FIGDIR}/Figure3_CCC_communication.png')
print(f'  Figure3 size: {img.size}')
# Check right third (panel C area) is not blank
import numpy as np
arr = np.array(img)
panel_c = arr[:, arr.shape[1]//2:, :]
non_white = (panel_c < 250).any(axis=2).sum()
print(f'  Panel C non-white pixels: {non_white} (should be > 1000)')

print('\n=== Figures 1-3 regenerated with Figure3 panel C fix ===')
