# Track 1: LASSO Biomarker Panel from OSDR Data

## Overview
Identifies a stable spaceflight biomarker gene panel from 6 NASA OSDR *Arabidopsis* transcriptomic studies using LASSO regression with leave-one-study-out cross-validation (LOOS-CV).

## Scripts
| Script | Description |
|--------|-------------|
| `04_generate_lasso_figures.py` | Generates Figure 1 (ROC/AUC, stability, coefficient forest) |
| `05_fix_figure3_panelC.py` | Fixes and regenerates Figure 3 panel C (top LR pairs) |

## Data
| File | Description |
|------|-------------|
| `biomarker_panel.csv` | 85 stable biomarker genes with coefficients |
| `all_feature_stability.csv` | All 2000 features with selection frequencies |
| `discovery_performance.csv` | Per-study AUC, sensitivity, specificity |
| `osdr_pooled_metadata.csv` | 156 sample metadata (flight/ground, study, tissue) |
| `osdr_study_catalog.csv` | 6 study descriptions |

## Key Results
- Mean CV AUC: 0.734 (95% CI: 0.431–1.0)
- 85 stable features at ≥50% selection frequency
- Top 4 genes (100% selection): AT4G10250, AT3G07365, AT2G14247, AT4G01390

## Input Data (external)
- OSDR studies: OSD-37, OSD-678, OSD-38, OSD-321, OSD-120, OSD-624
- Retrieved from: https://visualization.osdr.nasa.gov/biodata/api/
