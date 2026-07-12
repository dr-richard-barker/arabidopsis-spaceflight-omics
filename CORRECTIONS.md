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

## Pipeline sweep — completed 2026-07-11
- **Scripts fixed in place** (correct AGIs baked in): `05_ca2_pathway.R`, `08_ca2_k_pathway_network.R`
  (atomic single-pass AGI remap), `22_make_sankey_panel_b.py` (CBL1 title). The interim `*_corrected`
  copies were folded into the originals and removed.
- **Data tables regenerated** from the atlas with correct loci and distributed to all three locations
  (`track3_ccc_ggplantmap/data/`, `tables/`, `manuscript/supplementary/supplementary_tables/`):
  `ca2_k_expr_per_celltype.csv`, `ca2_k_coexpression_matrix.csv` + per-cell-type,
  `ca2_k_circuit_nodes.csv`, `ca2_k_signaling_circuit.csv` (via `13_circuit_summary_table.py`).
  `ca2_k_ccc_circuit.csv` is cell-type-level and was unaffected.
- **Reproducibility script added:** `track3_ccc_ggplantmap/scripts/29_regen_circuit_base_tables.py`
  regenerates the base expression/co-expression tables from the atlas AnnData.
- **Verified:** re-running the pipeline from the clean repo data reproduces the corrected figures
  (network = 27 nodes / 56 edges; CIPK23 guard-highest 0.333; guard co-expression module = CBL1–ANNAT1).
- The raw PlantCellChat database files (`ccc_lr_pairs.csv`, `ccc_full_communication_table.csv`) were
  **deliberately not modified** — they legitimately contain AT4G35310 (CPK5), AT1G30270 (CIPK23), etc.
  as distinct real entries; the error was only in the downstream alias mapping.

## Note
- Figure 8 is a raster (PNG) composite; the previous `Figure8_*.svg` (old ggKEGG version) was removed
  to avoid a stale vector. Vector sub-panels are available in `work/` if a vector Figure 8 is required.

## Unaffected
LASSO panel (Track 1), scPlantLLM (Track 2), and pathway-level CCC results (Ca²⁺ dominance 249.5,
10.2×, communication matrix) — these did not depend on the mislabeled loci.
