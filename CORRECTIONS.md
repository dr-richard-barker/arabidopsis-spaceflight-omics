# Corrections — gene-identity fix (Ca²⁺/K⁺ circuit)

**Date:** 2026-07-11

An audit found that the Ca²⁺/K⁺ circuit analysis used **incorrect AGI loci** for 14 of 23 named
genes (e.g. "CBL9" = AT5G24270 was actually SOS3/CBL4; "CIPK23" = AT4G35310 was actually CPK5;
the Shaker K⁺ channels were cyclically permuted). Correct loci are in
[`tables/corrected_gene_agi_map.csv`](tables/corrected_gene_agi_map.csv); full analysis in
`GENE_ID_AUDIT.md` (repository root of the working copy).

## What was corrected
- **Manuscript** (`manuscript/manuscript.docx`, `.pdf`, `manuscript_with_figures.docx`,
  `…_ACTIVE_DRAFT.docx`): correct AGIs; cascade reframed **CBL9 → CBL1** (the sensor actually present
  in the PlantCellChat network); the CIPK23–CaM3 / CBL9–AKT1 co-expression claims removed (they were
  CPK5 / SOS3 artifacts); the ggKEGG ath04075 panel removed (real CIPK23 is not in that pathway —
  the panel showed CPK5); network updated to 27 nodes / 56 edges.
- **Figures** re-generated from the atlas (GSE226097 seedling_6d, 41,314 cells) with correct loci:
  Figure 5, Figure 7 (all panels), the redesigned Figure 8 (leaf maps + Sankey, no KEGG), and
  Supplementary S4/S5/S8.
- **Scripts:** corrected re-run tooling added — `track3_ccc_ggplantmap/scripts/27_recompute_circuit_correct_ids.py`,
  `28_fig8_leaf_corrected.R`, and `*_corrected.R/.py` copies of the network/Sankey scripts (atomic AGI remap).

## Not fully cleaned yet (follow-up)
- The per-track **data tables** (`track3_ccc_ggplantmap/data/ca2_k_*.csv`) and the **original**
  scripts `05_ca2_pathway.R`, `08_ca2_k_pathway_network.R`, `22_make_sankey_panel_b.py` still contain
  the old AGIs; use the `*_corrected` versions + `corrected_gene_agi_map.csv` to regenerate. A full
  data/script sweep is recommended before final release.
- Figure 8 is a raster (PNG) composite; the previous `Figure8_*.svg` (old ggKEGG version) was removed
  to avoid a stale vector. Vector sub-panels are available in `work/` if a vector Figure 8 is required.

## Unaffected
LASSO panel (Track 1), scPlantLLM (Track 2), and pathway-level CCC results (Ca²⁺ dominance 249.5,
10.2×, communication matrix) — these did not depend on the mislabeled loci.
