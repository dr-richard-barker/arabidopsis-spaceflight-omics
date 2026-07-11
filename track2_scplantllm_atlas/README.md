# Track 2: scPlantLLM Zero-Shot Atlas Integration

## Overview
Applies the scPlantLLM plant foundation model in zero-shot mode to a 9,998-cell *Arabidopsis* seedling atlas, evaluating cell type annotation accuracy and clustering quality.

## Scripts
| Script | Description |
|--------|-------------|
| `01_build_anndata.py` | Converts Seurat RDS to AnnData h5ad format |
| `02_convert_atlas.R` | Atlas conversion (Seurat → AnnData) |
| `05_generate_figures.py` | Generates Figure 2 (UMAP panels) and Figure 4/5 assembly |

## Data
| File | Description |
|------|-------------|
| `scplantllm_zero_shot_annotations.csv` | Per-cell zero-shot predictions |
| `cluster_assignments.csv` | Leiden cluster labels |
| `cluster_degs.csv` | Top DEGs per cluster |
| `cluster_degs_all.csv` | All DEG results |
| `cluster_evaluation.csv` | ARI, NMI metrics |
| `cluster_vs_celltype_crosstab.csv` | Cluster × cell type contingency |
| `cluster_vs_zeroshot_crosstab.csv` | Cluster × zero-shot prediction contingency |
| `umap_coords.csv` | UMAP coordinates for 9,998 cells |

## Key Results
- Exact accuracy: 19.0%, Loose mapping accuracy: 28.1%
- Leiden res=0.5: 6 clusters, ARI=0.033, NMI=0.040

## Input Data (external)
- Seedling atlas: GEO GSE226097 (Lee et al., Nature Plants 2025)
- scPlantLLM pretrained weights: Cao et al., 2025
- Weight conversion: flash_attn Wqkv → PyTorch in_proj (see Supplementary Methods S4)
