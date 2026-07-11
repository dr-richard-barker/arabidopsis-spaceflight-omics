# Track 3: Cell-Cell Communication, Spatial Visualization, and GO Enrichment

## Overview
Characterizes intercellular signaling using PlantCellChat, visualizes expression on plant anatomy with ggPlantmap, constructs a Ca2+ pathway network, and performs GO enrichment analysis on the 85-gene LASSO panel.

## Scripts
| Script | Description |
|--------|-------------|
| `04_ggplantmap_composites.R` | ggPlantmap spatial visualization composites |
| `05_ca2_pathway.R` | Ca2+ pathway network (stress-majorization layout) |
| `06_gene_functions.py` | Gene function annotations for top 4 LASSO genes |
| `07_go_enrichment.R` | GO/KEGG enrichment via g:Profiler (custom background, BH-FDR) |
| `08_ca2_k_pathway_network.R` | Expanded Ca2+/K+ pathway network (27 nodes, 67 edges) with ggpathway |
| `09_ca2_k_coexpression.py` | Per-cell-type co-expression matrices for 25 Ca2+/K+ genes |
| `10_organ_marker_analysis.py` | Tissue-specific marker diagnostic for organ inference |
| `11_ca2_k_ccc_circuit.py` | Ca2+ and K+ cell-to-cell communication circuit |
| `12_ggplantmap_ca2_k_circuit.R` | ggPlantmap maps (seedling, root tip, leaf) for circuit genes |
| `13_circuit_summary_table.py` | Ca2+/K+ signaling circuit summary tables |
| `14_figure7_composite.py` | Figure 7 composite (network, heatmap, spatial map, CCC circuit) |
| `15_sankey_ca2_k_cascade.py` | Two-panel Sankey (cell-type flow + molecular cascade) |
| `16_ggplantmap_montage.py` | Supplementary ggPlantmap montage (3-panel stack) |
| `17_reassemble_figure5.py` | Figure 5 reassembly (refined 20-node network) |
| `18_ggkegg_figure8a.R` | **Figure 8 panel a**: CBL9-CIPK23-AKT1 cascade on KEGG ath04075 with guard cell expression overlay (ggkegg) |
| `19_combine_figure8.py` | Combine ggKEGG panel a + molecular cascade Sankey panel b into Figure 8 |

## Data
| File | Description |
|------|-------------|
| `ccc_communication_strength.csv` | 5×5 cell type communication strength matrix |
| `ccc_full_communication_table.csv` | All 72,350 LR interactions with probabilities |
| `ccc_lr_pairs.csv` | PlantCellChat LR database (3,140 interactions) |
| `ccc_lr_pairs_with_strength.csv` | 2,894 measured LR pairs with total communication |
| `ccc_per_signal_strength.csv` | Per-signal, per-cell-type communication strength |
| `ccc_signaling_strength.csv` | Total strength per signaling pathway |
| `ccc_diffexp_signaling.csv` | Differential signaling genes |
| `gene_function_annotations.csv` | Top 4 LASSO gene annotations (12 fields each) |
| `ca2_signaling_analysis.json` | Ca2+ pathway analysis summary |
| `go_enrichment_85_genes.csv` | Full GO enrichment results (1,481 terms) |
| `go_enrichment_significant.csv` | 14 significant terms (BH-FDR < 0.05) |
| `go_enrichment_summary.csv` | Summary with -log10(p) and gene lists |

## Key Results
- Ca2+ signaling dominates: strength = 249.5 (4.6× next pathway)
- GO enrichment: flight-up genes → heat/ER-proteostasis; ground-up → flavonoid metabolism
- Ca2+ network: 20 nodes, 43 edges, stress-majorization layout
- Ca2+/K+ crosstalk circuit: CBL9-CIPK23-AKT1 cascade, guard-cell-specific co-expression (CDPK4-ANNAT1 r = 0.56, CIPK23-CaM3 r = 0.46, CBL9-AKT1 r = 0.34)
- Ca2+ 10.2-fold stronger than K+ signaling; guard cells dominant Ca2+ source (outgoing strength 57.1)
- Figure 8 panel a maps the cascade onto KEGG ath04075 (Plant hormone signal transduction); genes absent from KEGG (CBL9, CBL1, CBL2, CML24, AKT1, KC1, ANNAT1) added as custom expression-colored nodes

## Input Data (external)
- PlantCellChat: https://github.com/mrliuw/PlantCellChat
- ggPlantmap: https://github.com/gojo-jojo/ggPlantmap (Jo & Kajala, J Exp Bot 2024)
- ggkegg: https://github.com/noriakis/ggkegg (KEGG pathway ath04075)
- g:Profiler: https://biit.cs.ut.ee/gprofiler (organism: athaliana)
