# Arabidopsis Spaceflight Omics: Biomarker Discovery and Single-Cell Atlas Integration

A multi-track computational pipeline for identifying spaceflight-responsive biomarker genes in *Arabidopsis thaliana* using NASA OSDR transcriptomic data, validating them against a single-cell seedling atlas, and characterizing cell-cell communication patterns.

## Overview

This repository contains the complete analysis pipeline for a study submitted to npj Microgravity. The work integrates three complementary approaches:

1. **Track 1 — LASSO Biomarker Panel (OSDR)**: Identifies a stable 85-gene biomarker panel predictive of spaceflight condition from 6 NASA OSDR studies (156 samples) using LASSO regression with leave-one-study-out cross-validation.

2. **Track 2 — scPlantLLM Atlas Integration**: Applies the scPlantLLM plant foundation model in zero-shot mode to a 9,998-cell *Arabidopsis* seedling atlas, evaluating cell type annotation accuracy and clustering quality.

3. **Track 3 — Cell-Cell Communication & Spatial Visualization**: Characterizes intercellular signaling using PlantCellChat, with focus on Ca2+ pathway dominance, and visualizes expression patterns on plant anatomy using ggPlantmap.

## Repository Structure

```
arabidopsis-spaceflight-omics/
├── track1_lasso_osdr/          # LASSO biomarker discovery from OSDR data
│   ├── scripts/                # Python/R analysis scripts
│   ├── data/                   # Biomarker panel, stability, performance tables
│   └── figures/                # ROC, stability, coefficient figures
├── track2_scplantllm_atlas/    # scPlantLLM zero-shot + clustering
│   ├── scripts/                # Atlas conversion, weight conversion, inference
│   ├── data/                   # Annotations, cluster assignments, DEGs
│   └── figures/                # UMAP, clustering metrics, DEG heatmap
├── track3_ccc_ggplantmap/      # Cell-cell communication & spatial viz
│   ├── scripts/                # PlantCellChat, ggPlantmap, Ca2+ pathway, GO enrichment
│   ├── data/                   # CCC tables, gene annotations, GO enrichment
│   └── figures/                # CCC heatmaps, ggPlantmap composites, Ca2+ network
├── manuscript/                 # Manuscript and supplementary materials
│   ├── manuscript.docx
│   ├── manuscript.pdf
│   ├── Supplementary_Methods.md
│   ├── figures/                # 8 main figures (PNG + SVG)
│   └── supplementary/          # Supplementary figures and tables
├── integration/                # Cross-track integration analysis
│   ├── integration_analysis.py
│   └── integration_summary.json
├── requirements.txt            # Python dependencies
├── environment.R               # R dependencies
├── LICENSE                     # MIT License
├── CITATION.cff                # Zenodo citation metadata
└── .gitignore
```

## Key Results

- **LASSO panel**: Mean CV AUC = 0.734, 85 stable biomarker genes, top 4 selected at 100% frequency
- **scPlantLLM**: Zero-shot loose accuracy = 28.1%, highlighting domain shift challenges for plant foundation models
- **Cell-cell communication**: Ca2+ signaling dominates (strength = 249.5, 4.6× next pathway)
- **GO enrichment**: Flight-up genes enriched for heat/ER-proteostasis; ground-up genes enriched for flavonoid metabolism
- **Ca2+/K+ crosstalk circuit**: CBL9-CIPK23-AKT1 cascade connects Ca2+ sensing to K+ uptake; guard cells are the dominant Ca2+ source (10.2-fold over K+); strongest guard cell co-expression: CDPK4-ANNAT1 (r=0.56), CBL9-AKT1 (r=0.34)
- **Sankey diagram**: Two-panel visualization tracing signal flow at cell-type and molecular cascade levels
- **Spatial mapping**: ggPlantmap montage showing circuit gene expression across seedling, root tip, and leaf anatomy

## Data Sources

- **NASA OSDR**: https://visualization.osdr.nasa.gov/biodata/api/ (studies OSD-37, OSD-678, OSD-38, OSD-321, OSD-120, OSD-624)
- **Seedling atlas**: GEO GSE226097 (Lee et al., Nature Plants 2025)
- **scPlantLLM**: Cao et al., Genomics Proteomics Bioinformatics 2025

## Installation

### Python dependencies
```bash
pip install -r requirements.txt
```

### R dependencies
```r
# Install from CRAN
install.packages(c("igraph", "ggraph", "ggplot2", "dplyr", "tidyr", "viridis",
                   "cowplot", "patchwork", "gprofiler2", "remotes"))
# Install from Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "GO.db"))
# ggPlantmap from GitHub
remotes::install_github("leonardojo/ggPlantmap")
# PlantCellChat from GitHub
remotes::install_github("mrliuw/PlantCellChat")
```

## Reproducibility & paths

All scripts resolve their input/output locations **relative to the repository root** (computed from each script's own location), so they run from any working directory with no edits. The layout each script expects:

| Base | Default location | Holds |
|------|------------------|-------|
| `RESULTS` | repo root (contains `tables/`, `figures/`) | input tables + generated tables/figures |
| `ATLAS` | `<repo>/atlas/` | large single-cell intermediates (not shipped — regenerate via `track2_scplantllm_atlas/scripts/01_build_anndata.py`) |
| `WORK` | `<repo>/work/` | scratch/intermediate outputs (auto-created) |

Override any base with an environment variable if your data lives elsewhere:

```bash
export ASO_ROOT=/path/to/results      # dir containing tables/ and figures/
export ASO_ATLAS=/path/to/atlas       # large .h5ad/.rds intermediates
export ASO_WORK=/path/to/scratch
```

The `tables/` directory (inputs) is version-controlled. The root `figures/` directory is git-ignored and **regenerated** by the analysis scripts — an empty `figures/` (with a `.gitkeep`) ships so scripts can write to it on a fresh clone. `atlas/`, `work/`, and `results/` are also git-ignored. (The publication figures under `manuscript/figures/` and each track's `figures/` remain version-controlled.)

## Usage

Each track can be run independently. See the README.md in each track directory for specific instructions.

## License

MIT License — see LICENSE file for details.

## Citation

See CITATION.cff for citation information.
