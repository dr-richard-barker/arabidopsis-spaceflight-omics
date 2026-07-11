#!/usr/bin/env python3
"""Build AnnData from HDF5 sparse matrix + CSV metadata"""
# --- portable paths (de-sandboxed; replaces /mnt/results and /workspace) ---
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
REPO_ROOT = _os.path.abspath(_os.path.join(_HERE, '..', '..'))
RESULTS = _os.environ.get("ASO_ROOT", REPO_ROOT)          # holds tables/ and figures/
ATLAS   = _os.environ.get("ASO_ATLAS", _os.path.join(REPO_ROOT, "atlas"))  # large intermediates (not shipped)
WORK    = _os.environ.get("ASO_WORK", _os.path.join(REPO_ROOT, "work"))    # scratch outputs
_os.makedirs(WORK, exist_ok=True)
# --- end portable paths ---

import h5py
import numpy as np
import scipy.sparse as sp
import pandas as pd
import anndata as ad
import scanpy as sc

print("=== Building AnnData from atlas HDF5 ===")

# Read sparse matrix from HDF5
print("  Reading sparse count matrix...")
f = h5py.File(ATLAS + '/seedling_6d.h5', 'r')
data = f['X/data'][:]
indices = f['X/indices'][:]
indptr = f['X/indptr'][:]
shape = tuple(f['X/shape'][:])
f.close()

counts = sp.csr_matrix((data, indices, indptr), shape=shape)
print(f"  Count matrix: {counts.shape[0]} genes x {counts.shape[1]} cells")
print(f"  Non-zero entries: {counts.nnz}")
print(f"  Density: {counts.nnz / (counts.shape[0] * counts.shape[1]):.4f}")

# Read gene names
genes = pd.read_csv(ATLAS + '/gene_names.csv')['gene'].values
print(f"  Genes: {len(genes)}")

# Read cell metadata
meta = pd.read_csv(ATLAS + '/cell_metadata.csv', index_col='cell_id')
print(f"  Cell metadata: {meta.shape}")
print(f"  Cell types: {meta['CellType'].value_counts().to_dict()}")

# Read UMAP
umap = pd.read_csv(ATLAS + '/umap_coords.csv', index_col='cell_id')
print(f"  UMAP coords: {umap.shape}")

# Build AnnData
# Note: counts matrix is genes x cells, AnnData expects cells x genes
adata = ad.AnnData(
    X=counts.T.tocsr(),  # transpose to cells x genes
    obs=meta,
    var=pd.DataFrame(index=genes),
)
adata.obsm['X_umap'] = umap.loc[adata.obs_names, ['UMAP_1', 'UMAP_2']].values

print(f"\n  AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  obs columns: {list(adata.obs.columns)}")

# QC
print("\n=== QC ===")
sc.pp.filter_cells(adata, min_genes=200)
print(f"  After min_genes=200 filter: {adata.shape[0]} cells")
sc.pp.filter_genes(adata, min_cells=3)
print(f"  After min_cells=3 filter: {adata.shape[1]} genes")

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3', layer=None, subset=False)
print(f"  HVGs: {adata.var['highly_variable'].sum()}")

# Save
adata.write_h5ad(ATLAS + '/seedling_6d_anndata.h5ad')
print(f"\n✓ Saved AnnData to " + ATLAS + "/seedling_6d_anndata.h5ad")
print(f"  Final: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  Cell types: {adata.obs['CellType'].value_counts().to_dict()}")
