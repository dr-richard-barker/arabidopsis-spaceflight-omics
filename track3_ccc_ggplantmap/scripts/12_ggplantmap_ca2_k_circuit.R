#!/usr/bin/env Rscript
# ggPlantmap spatial visualization of Ca2+/K+ circuit
# Maps gene expression onto: (1) seedling map, (2) root tip longitudinal, (3) leaf cross-section
# --- portable paths (de-sandboxed; replaces /mnt/results and /workspace) ---
.aso_here <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))
  getwd()
})
REPO_ROOT <- normalizePath(file.path(.aso_here, "..", ".."), mustWork = FALSE)
RESULTS <- Sys.getenv("ASO_ROOT",  unset = REPO_ROOT)          # holds tables/ and figures/
ATLAS   <- Sys.getenv("ASO_ATLAS", unset = file.path(REPO_ROOT, "atlas"))  # large intermediates (not shipped)
WORK    <- Sys.getenv("ASO_WORK",  unset = file.path(REPO_ROOT, "work"))    # scratch outputs
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
# --- end portable paths ---

invisible(NULL)  # removed sandbox-specific .Rlib path (use environment.R)
library(ggPlantmap)
library(ggplot2)
library(dplyr)
library(cowplot)
library(svglite)

cat("=== ggPlantmap Ca2+/K+ Circuit Visualization ===\n\n")

# Load expression data
expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"), stringsAsFactors = FALSE)

# Cell type -> organ/ROI mapping
# Seedling map: Cotyledon, Hook, Hypocotyl, Root (control only)
# Mesophyll+Guard -> Cotyledon, Epidermal -> Hypocotyl, Stele -> Root, Meristematic -> Root
ct_to_seedling <- data.frame(
  ROI.name = c("Cotyledon_control", "Hook_control", "Hypocotyl_control", "Root_control"),
  Epidermal = c(0, 0, 0, 0),  # Epidermal is shared, use average
  Mesophyll = c(0, 0, 0, 0),
  Stele = c(0, 0, 0, 0),
  Meristematic = c(0, 0, 0, 0),
  Guard = c(0, 0, 0, 0)
)
# Map: Cotyledon = Mesophyll + Guard avg; Hypocotyl = Epidermal + Stele avg; Root = Stele + Meristematic avg; Hook = avg of all
for (gene_alias in unique(expr$alias)) {
  gene_expr <- expr %>% filter(alias == gene_alias)
  cotyledon_val <- mean(c(gene_expr$mean_expr[gene_expr$celltype == "Mesophyll"],
                          gene_expr$mean_expr[gene_expr$celltype == "Guard"]), na.rm = TRUE)
  hypocotyl_val <- mean(c(gene_expr$mean_expr[gene_expr$celltype == "Epidermal"],
                          gene_expr$mean_expr[gene_expr$celltype == "Stele"]), na.rm = TRUE)
  root_val <- mean(c(gene_expr$mean_expr[gene_expr$celltype == "Stele"],
                     gene_expr$mean_expr[gene_expr$celltype == "Meristematic"]), na.rm = TRUE)
  hook_val <- mean(c(cotyledon_val, hypocotyl_val, root_val), na.rm = TRUE)
  ct_to_seedling[ct_to_seedling$ROI.name == "Cotyledon_control", as.character(gene_alias)] <- cotyledon_val
  ct_to_seedling[ct_to_seedling$ROI.name == "Hypocotyl_control", as.character(gene_alias)] <- hypocotyl_val
  ct_to_seedling[ct_to_seedling$ROI.name == "Root_control", as.character(gene_alias)] <- root_val
  ct_to_seedling[ct_to_seedling$ROI.name == "Hook_control", as.character(gene_alias)] <- hook_val
}

# Key circuit genes to visualize
circuit_genes <- c("AKT1", "CBL9", "CIPK23", "CML24", "KC1", "CDPK3")

# ============================================================
# 1. Seedling map: circuit genes across cotyledon/hypocotyl/root
# ============================================================
cat("=== Seedling map ===\n")
seedling_map <- ggPm.At.seedling.saltdrought %>% filter(treatment == "control") %>% group_by(ROI.name, ROI.id, organ, treatment) %>% mutate(point = row_number()) %>% ungroup()

plots_seedling <- list()
for (gene in circuit_genes) {
  values <- ct_to_seedling %>% select(ROI.name, all_of(gene))
  colnames(values)[2] <- "value"
  merged <- ggPlantmap.merge(seedling_map, values, "ROI.name")
  p <- ggPlantmap.heatmap(merged, value.quant = value) +
    scale_fill_gradient(low = "white", high = "darkred", name = "Expr") +
    ggtitle(gene) +
    theme_void() + coord_fixed() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  plots_seedling[[gene]] <- p
}

# Composite: 6 genes in 2x3 grid
seedling_composite <- plot_grid(plotlist = plots_seedling, ncol = 3, labels = c("a","b","c","d","e","f"),
                                label_size = 12, align = "hv")
title <- ggdraw() + draw_label("Ca2+/K+ Circuit Genes on Seedling Anatomy (cell type proxy)",
                               fontface = "bold", size = 13)
seedling_final <- plot_grid(title, seedling_composite, ncol = 1, rel_heights = c(0.1, 1))

ggsave(paste0(RESULTS, "/figures/ggplantmap_ca2_k_circuit_seedling.png"), seedling_final,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_ca2_k_circuit_seedling.svg"), seedling_final,
       width = 14, height = 10, bg = "white")
cat("Saved: ggplantmap_ca2_k_circuit_seedling.png/.svg\n")

# ============================================================
# 2. Root tip longitudinal: Ca2+ cascade in root zones
# ============================================================
cat("=== Root tip longitudinal ===\n")
root_map <- ggPm.At.roottip.longitudinal

# Map cell types to root zones
# Meristematic -> Meristem zones, Stele -> Stele zones, Epidermal -> Epidermis
root_zone_expr <- data.frame(
  ROI.name = unique(root_map$ROI.name),
  AKT1 = 0, CBL9 = 0, CIPK23 = 0, CML24 = 0, KC1 = 0, CDPK3 = 0
)
for (roi in root_zone_expr$ROI.name) {
  if (grepl("Meristem", roi) || grepl("Columella", roi) || grepl("Lateral", roi)) {
    ct <- "Meristematic"
  } else if (grepl("Stele", roi) || grepl("Vascular", roi) || grepl("Pericycle", roi)) {
    ct <- "Stele"
  } else if (grepl("Epidermis", roi)) {
    ct <- "Epidermal"
  } else if (grepl("Cortex", roi) || grepl("Endodermis", roi)) {
    ct <- "Mesophyll"  # closest proxy for cortex
  } else {
    ct <- "Meristematic"
  }
  for (gene in circuit_genes) {
    val <- expr$mean_expr[expr$alias == gene & expr$celltype == ct]
    if (length(val) > 0) root_zone_expr[root_zone_expr$ROI.name == roi, gene] <- val
  }
}

plots_root <- list()
for (gene in circuit_genes) {
  values <- root_zone_expr %>% select(ROI.name, all_of(gene))
  colnames(values)[2] <- "value"
  merged <- ggPlantmap.merge(root_map, values, "ROI.name")
  p <- ggPlantmap.heatmap(merged, value.quant = value) +
    scale_fill_gradient(low = "white", high = "darkred", name = "Expr") +
    ggtitle(gene) +
    theme_void() + coord_fixed() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  plots_root[[gene]] <- p
}

root_composite <- plot_grid(plotlist = plots_root, ncol = 3, labels = c("a","b","c","d","e","f"),
                            label_size = 12, align = "hv")
title2 <- ggdraw() + draw_label("Ca2+/K+ Circuit Genes on Root Tip Longitudinal Section",
                                fontface = "bold", size = 13)
root_final <- plot_grid(title2, root_composite, ncol = 1, rel_heights = c(0.1, 1))

ggsave(paste0(RESULTS, "/figures/ggplantmap_ca2_k_circuit_roottip.png"), root_final,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_ca2_k_circuit_roottip.svg"), root_final,
       width = 14, height = 10, bg = "white")
cat("Saved: ggplantmap_ca2_k_circuit_roottip.png/.svg\n")

# ============================================================
# 3. Leaf cross-section: K+ channels in leaf tissues
# ============================================================
cat("=== Leaf cross-section ===\n")
leaf_map <- ggPm.At.leaf.crosssection

# Map cell types to leaf ROIs
leaf_roi_expr <- data.frame(
  ROI.name = unique(leaf_map$ROI.name),
  AKT1 = 0, CBL9 = 0, CIPK23 = 0, CML24 = 0, KC1 = 0, CDPK3 = 0
)
for (roi in leaf_roi_expr$ROI.name) {
  if (grepl("stomata", roi)) {
    ct <- "Guard"
  } else if (grepl("epidermis", roi)) {
    ct <- "Epidermal"
  } else if (grepl("Parenchima", roi)) {
    ct <- "Mesophyll"
  } else if (grepl("vascular", roi)) {
    ct <- "Stele"
  } else {
    ct <- "Mesophyll"
  }
  for (gene in circuit_genes) {
    val <- expr$mean_expr[expr$alias == gene & expr$celltype == ct]
    if (length(val) > 0) leaf_roi_expr[leaf_roi_expr$ROI.name == roi, gene] <- val
  }
}

plots_leaf <- list()
for (gene in circuit_genes) {
  values <- leaf_roi_expr %>% select(ROI.name, all_of(gene))
  colnames(values)[2] <- "value"
  merged <- ggPlantmap.merge(leaf_map, values, "ROI.name")
  p <- ggPlantmap.heatmap(merged, value.quant = value) +
    scale_fill_gradient(low = "white", high = "darkred", name = "Expr") +
    ggtitle(gene) +
    theme_void() + coord_fixed() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  plots_leaf[[gene]] <- p
}

leaf_composite <- plot_grid(plotlist = plots_leaf, ncol = 3, labels = c("a","b","c","d","e","f"),
                            label_size = 12, align = "hv")
title3 <- ggdraw() + draw_label("Ca2+/K+ Circuit Genes on Leaf Cross-Section",
                                fontface = "bold", size = 13)
leaf_final <- plot_grid(title3, leaf_composite, ncol = 1, rel_heights = c(0.1, 1))

ggsave(paste0(RESULTS, "/figures/ggplantmap_ca2_k_circuit_leaf.png"), leaf_final,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_ca2_k_circuit_leaf.svg"), leaf_final,
       width = 14, height = 10, bg = "white")
cat("Saved: ggplantmap_ca2_k_circuit_leaf.png/.svg\n")

cat("\nAll ggPlantmap circuit visualizations complete.\n")
