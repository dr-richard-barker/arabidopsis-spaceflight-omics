#!/usr/bin/env python3
"""
Re-compute the Ca2+/K+ circuit per-cell-type expression and co-expression using
the CORRECT gene->AGI identifiers (see GENE_ID_AUDIT.md / corrected_gene_agi_map.csv).

The original circuit tables were built from a scrambled name->AGI alias table
(14 of 23 genes wrong; e.g. "CBL9"=AT5G24270 was actually SOS3/CBL4,
"CIPK23"=AT4G35310 was actually CPK5). This script re-derives the expression and
co-expression matrices from the atlas for the correct loci.

REQUIRES the atlas AnnData (not shipped): ATLAS/seedling_6d_anndata.h5ad
(cells x genes, log-normalized; var_names = AGI; obs['CellType']).
Build it first with track2_scplantllm_atlas/scripts/01_build_anndata.py, or set
ASO_ATLAS to a directory containing seedling_6d_anndata.h5ad.

Outputs (written to WORK so originals are not overwritten until reviewed):
  WORK/ca2_k_expr_per_celltype_CORRECTED.csv
  WORK/ca2_k_coexpression_<CellType>_CORRECTED.csv  (one per cell type)
  WORK/coexpr_before_after.csv                       (headline pairs, old vs new)
"""
# --- portable paths ---
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
REPO_ROOT = _os.path.abspath(_os.path.join(_HERE, '..', '..'))
RESULTS = _os.environ.get("ASO_ROOT", REPO_ROOT)
ATLAS   = _os.environ.get("ASO_ATLAS", _os.path.join(REPO_ROOT, "atlas"))
WORK    = _os.environ.get("ASO_WORK", _os.path.join(REPO_ROOT, "work"))
_os.makedirs(WORK, exist_ok=True)
# --- end portable paths ---

import sys
import numpy as np
import pandas as pd

MAP_CSV = _os.environ.get("ASO_GENEMAP",
                          _os.path.join(REPO_ROOT, "tables", "corrected_gene_agi_map.csv"))
ATLAS_H5AD = _os.path.join(ATLAS, "seedling_6d_anndata.h5ad")

def die(msg):
    print("ERROR:", msg); sys.exit(1)

if not _os.path.exists(MAP_CSV):
    die(f"corrected gene->AGI map not found: {MAP_CSV}")
if not _os.path.exists(ATLAS_H5AD):
    die(f"atlas AnnData not found: {ATLAS_H5AD}\n"
        f"       Build it with track2 01_build_anndata.py or set ASO_ATLAS.")

gmap = pd.read_csv(MAP_CSV)
# name -> correct AGI
name2agi = dict(zip(gmap['gene'], gmap['correct_AGI']))
print(f"Loaded {len(name2agi)} corrected gene->AGI mappings")

import anndata as ad
adata = ad.read_h5ad(ATLAS_H5AD)
print(f"Atlas: {adata.shape[0]} cells x {adata.shape[1]} genes; "
      f"cell types = {sorted(adata.obs['CellType'].unique())}")

var = set(map(str, adata.var_names))
present = {n: a for n, a in name2agi.items() if a in var}
missing = {n: a for n, a in name2agi.items() if a not in var}
if missing:
    print(f"\nWARNING: {len(missing)} correct AGIs absent from the atlas "
          f"(cannot be quantified): "
          + ", ".join(f"{n}({a})" for n, a in missing.items()))
genes = list(present.keys())
agis = [present[n] for n in genes]

# dense expression sub-matrix (cells x selected genes)
X = adata[:, agis].X
X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
expr_df = pd.DataFrame(X, columns=genes, index=adata.obs_names)
expr_df['CellType'] = adata.obs['CellType'].values

# --- per-cell-type mean expression + pct expressing ---
rows = []
for ct, sub in expr_df.groupby('CellType'):
    m = sub[genes]
    for g in genes:
        rows.append({'gene': present[g], 'alias': g, 'celltype': ct,
                     'mean_expr': float(m[g].mean()),
                     'pct_expr': float((m[g] > 0).mean() * 100)})
    # ALL row
for g in genes:
    rows.append({'gene': present[g], 'alias': g, 'celltype': 'ALL',
                 'mean_expr': float(expr_df[g].mean()),
                 'pct_expr': float((expr_df[g] > 0).mean() * 100)})
expr_out = pd.DataFrame(rows)
expr_out.to_csv(_os.path.join(WORK, "ca2_k_expr_per_celltype_CORRECTED.csv"), index=False)
print(f"\nWrote corrected per-cell-type expression ({len(expr_out)} rows)")

# --- per-cell-type Pearson co-expression matrices ---
headline = [('CDPK4', 'ANNAT1'), ('CIPK23', 'CaM3'), ('CBL1', 'AKT1'),
            ('CBL9', 'AKT1'), ('CaM3', 'CDPK4'), ('CML24', 'KC1')]
ba = []
for ct, sub in expr_df.groupby('CellType'):
    corr = sub[genes].corr(method='pearson')
    corr.to_csv(_os.path.join(WORK, f"ca2_k_coexpression_{ct}_CORRECTED.csv"))
    if ct == 'Guard':
        for a, b in headline:
            v = corr.loc[a, b] if (a in corr.index and b in corr.columns) else np.nan
            ba.append({'cell_type': ct, 'pair': f"{a}-{b}",
                       'r_corrected_genes': round(float(v), 3) if v == v else 'n/a'})
pd.DataFrame(ba).to_csv(_os.path.join(WORK, "coexpr_before_after.csv"), index=False)
print("Wrote per-cell-type corrected co-expression matrices")
print("\nGuard-cell co-expression for the corrected loci (compare to the "
      "manuscript's current values):")
for r in ba:
    print(f"  {r['pair']:16} r = {r['r_corrected_genes']}")
print(f"\nAll outputs in: {WORK}")
print("Review, then (if sound) copy the *_CORRECTED tables over the originals and "
      "re-run scripts 09-26 to regenerate Figures 5/7/8 + supplements.")
