# npj Microgravity Manuscript Package

## Title
LASSO biomarker panel and single-cell atlas integration reveals Ca2+-dominated cell-cell communication in Arabidopsis spaceflight response

## Contents

### Main Manuscript
- `manuscript.docx` — Word format manuscript (all sections; ~5,700 words incl. figure legends)
- `manuscript.pdf` — PDF version of the manuscript
- `manuscript_with_figures.docx` — Word manuscript with all 8 figures embedded inline (review-friendly)
- `manuscript_with_figures_ACTIVE_DRAFT.docx` — **active working draft** (in progress, ~6,900 words); expanded text under revision. The frozen `manuscript.docx` / `manuscript_with_figures.docx` remain the reference submission versions until this draft is finalised.

### Main Figures (`figures/`)
8 main figures in PNG (300 dpi) and SVG (editable vector) formats:

1. **Figure1_LASSO_panel** — LASSO biomarker panel: (a) ROC/AUC per study, (b) feature stability, (c) coefficient forest plot
2. **Figure2_scPlantLLM_clustering** — UMAP: (a) atlas cell types, (b) Leiden clusters, (c) zero-shot predictions
3. **Figure3_CCC_communication** — Cell-cell communication: (a) heatmap, (b) top signaling pathways, (c) top LR pairs
4. **Figure4_ggPlantmap_LASSO_genes** — Spatial expression of 4 LASSO genes on (a) leaf and (b) root tip cross-sections
5. **Figure5_Ca2_pathway** — Ca2+ signaling: (a) data-driven ggpathway network, (b) conceptual schematic
6. **Figure6_integration** — Integration: (a) LASSO panel score on leaf, (b) Ca2+ signaling on leaf, (c) gene function table
7. **Figure7_ca2_k_circuit_composite** — Ca2+/K+ crosstalk signaling circuit
8. **Figure8_sankey_ca2_k_cascade** — CBL9–CIPK23–AKT1 cascade: KEGG systems biology context (ggKEGG overlay on ath04075), spatial expression, and molecular cascade (Sankey)

### Supplementary Materials (`supplementary/`)
- **`supplementary_figures/`** — 48 supplementary figure files (PNG + SVG) including:
  - Original LASSO plots (ROC curve, AUC distribution, stability barplot, coefficient forest, feature heatmap, calibration curve)
  - Original scPlantLLM plots (UMAP celltype, UMAP leiden, UMAP zeroshot, DEG heatmap, clustering metrics)
  - Original CCC plots (chord diagram, heatmap, LR pair stats, top LR pairs, top signals)
  - ggPlantmap composites (5 composites: LASSO genes on leaf/root, CCC communication, Ca2+ signaling, panel score)
  - Individual ggPlantmap maps (8 original maps)
  - Ca2+ pathway diagrams (ggpathway network, conceptual schematic)

- **`supplementary_tables/`** — 26 supplementary data tables (CSV + JSON) including:
  - `biomarker_panel.csv` — 85 stable LASSO features with coefficients
  - `all_feature_stability.csv` — All 2000 features ranked by selection frequency
  - `discovery_performance.csv` — Per-fold LASSO AUC/sensitivity/specificity
  - `gene_function_annotations.csv` — Top 4 LASSO gene function annotations
  - `ca2_signaling_analysis.json` — Ca2+ signaling pathway analysis
  - `ccc_communication_strength.csv` — 5×5 communication matrix
  - `ccc_full_communication_table.csv` — 72,350-row full communication table
  - `ccc_per_signal_strength.csv` — Per-signal per-cell-pair strength
  - `cluster_assignments.csv` — 9,998 cell cluster assignments
  - `cluster_degs.csv` — 260 filtered DEGs
  - `cluster_evaluation.csv` — ARI/NMI metrics
  - `scplantllm_zero_shot_annotations.csv` — 9,998 cell predictions
  - `umap_coords.csv` — UMAP coordinates
  - `integration_summary.json` — Complete integration summary
  - `osdr_pooled_metadata.csv` — 156 sample metadata
  - `osdr_study_catalog.csv` — Study catalog

## Key Results Summary

| Metric | Value |
|--------|-------|
| LASSO mean cross-study AUC | 0.734 |
| Stable features (≥50% selection) | 85 |
| Top 4 genes (100% selection) | AT4G10250, AT3G07365, AT2G14247, AT4G01390 |
| Atlas cells analyzed | 9,998 (subsampled from 41,314) |
| Zero-shot accuracy (loose mapping) | 28.1% |
| Leiden ARI vs cell type | 0.033 |
| Dominant CCC pathway | Ca2+ (strength 249.5, 4.6× next pathway) |
| CCC LR interactions | 3,140 |
| LASSO genes in atlas | 74/85 |
| LASSO genes as CCC ligands | 3 |

## Software and Data Sources
- NASA OSDR: visualization.osdr.nasa.gov/biodata/api/
- Seed-to-seed atlas: GEO GSE226097
- scPlantLLM: github.com/compbioNJU/scPlantLLM
- PlantCellChat: github.com/mrliuw/PlantCellChat
- ggPlantmap: github.com/leonardojo/ggPlantmap
- ggpathway: github.com/cxli233/ggpathway

## Notes
- Manuscript formatted for npj Microgravity initial submission (PDF or Word accepted)
- Figures at 300 dpi, Liberation Sans/Arial font, colorblind-friendly palettes
- Figure legends <350 words each, placed after references per npj Microgravity guidelines
- Results flagged as discovery-only (no external validation cohort)
- AT3G07365 is a natural antisense transcript (non-coding RNA), not a protein-coding gene
