#!/usr/bin/env Rscript
# ggPlantmap v5 - Redesigned using ggPlantmap.merge() + ggPlantmap.heatmap()
# Following the official ggPlantmap practical guide for single-cell data
# Creates polished composite figures matching the example repo aesthetics

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
library(tidyr)
library(viridis)
library(svglite)

cat("=== ggPlantmap v5 - Redesigned Composites ===\n\n")

# Load data
biomarker <- read.csv(paste0(RESULTS, "/tables/biomarker_panel.csv"))
ccc <- read.csv(paste0(RESULTS, "/tables/ccc_communication_strength.csv"), row.names = 1, check.names = FALSE)
signal_comm <- read.csv(paste0(RESULTS, "/tables/ccc_per_signal_strength.csv"), stringsAsFactors = FALSE)
avg_expr <- read.csv(paste0(RESULTS, "/tables/lasso_avg_expr_per_celltype.csv"), row.names = 1, check.names = FALSE)

# Top 4 LASSO genes
top_genes <- c("AT4G10250", "AT3G07365", "AT2G14247", "AT4G01390")

# Cell type -> ROI mapping
# Leaf cross-section
ct_to_roi_leaf <- list(
  "Epidermal" = c("epidermis.adaxial", "epidermis.abaxial"),
  "Mesophyll" = c("Parenchima.palisade", "Parenchima.sponge"),
  "Stele" = c("vascularbundle.xylem", "vascularbundle.phloem", "vascularbundle.bundlesheet"),
  "Guard" = c("epidermis.stomata"),
  "Meristematic" = c()  # no direct leaf ROI
)

# Root tip cross-section
ct_to_roi_root <- list(
  "Epidermal" = c("Epidermis"),
  "Mesophyll" = c("Cortex"),  # tissue analogy
  "Stele" = c("Procambium", "Xylem", "Phloem", "Pericycle"),
  "Meristematic" = c("Columella"),  # tissue analogy
  "Guard" = c()  # no root ROI
)

# Helper: build ROI -> value mapping from cell-type data
build_roi_values <- function(celltype_values, ct_to_roi) {
  roi_vals <- data.frame(ROI.name = character(), value = numeric(), stringsAsFactors = FALSE)
  for (ct in names(celltype_values)) {
    rois <- ct_to_roi[[ct]]
    if (is.null(rois) || length(rois) == 0) next
    for (roi in rois) {
      roi_vals <- rbind(roi_vals, data.frame(ROI.name = roi, value = celltype_values[ct], stringsAsFactors = FALSE))
    }
  }
  return(roi_vals)
}

# Helper: create polished ggPlantmap heatmap
make_heatmap <- function(map_data, roi_values, title, fill_name = "Expression",
                         fill_scale = "gradient", limits = NULL,
                         legend_position = "right") {
  # Merge using ggPlantmap.merge
  merged <- ggPlantmap.merge(map_data, roi_values, id.x = "ROI.name", id.y = "ROI.name")

  # Build plot using ggPlantmap.heatmap
  p <- ggPlantmap.heatmap(merged, value.quant = value) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(family = "sans", size = 12, face = "bold", hjust = 0.5),
      legend.title = element_text(family = "sans", size = 10, face = "bold"),
      legend.text = element_text(family = "sans", size = 9),
      legend.position = legend_position,
      plot.margin = margin(5, 5, 5, 5)
    )

  if (fill_scale == "gradient") {
    if (!is.null(limits)) {
      p <- p + scale_fill_gradient(low = "white", high = "#D32F2F",
                                    name = fill_name, limits = limits, na.value = "grey90")
    } else {
      p <- p + scale_fill_gradient(low = "white", high = "#D32F2F",
                                    name = fill_name, na.value = "grey90")
    }
  } else if (fill_scale == "viridis") {
    if (!is.null(limits)) {
      p <- p + scale_fill_viridis_c(name = fill_name, limits = limits, na.value = "grey90", option = "C")
    } else {
      p <- p + scale_fill_viridis_c(name = fill_name, na.value = "grey90", option = "C")
    }
  } else if (fill_scale == "diverging") {
    p <- p + scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#FF9400",
                                  midpoint = 0, name = fill_name, na.value = "grey90")
  }

  return(p)
}

# Load maps
data("ggPm.At.leaf.crosssection", package = "ggPlantmap")
data("ggPm.At.roottip.crosssection", package = "ggPlantmap")
leaf_map <- ggPm.At.leaf.crosssection
root_map <- ggPm.At.roottip.crosssection

cat("Leaf ROIs:", paste(unique(leaf_map$ROI.name), collapse = ", "), "\n")
cat("Root ROIs:", paste(unique(root_map$ROI.name), collapse = ", "), "\n\n")

# ============================================================
# COMPOSITE 1: 4 LASSO genes on leaf cross-section (faceted)
# ============================================================
cat("=== Composite 1: 4 LASSO genes on leaf ===\n")

# Build expression data for each gene
gene_plots_leaf <- list()
for (gene in top_genes) {
  if (gene %in% rownames(avg_expr)) {
    ct_vals <- as.numeric(avg_expr[gene, ])
    names(ct_vals) <- colnames(avg_expr)
  } else {
    # Gene not in atlas (e.g., AT3G07365 is non-coding)
    ct_vals <- setNames(rep(NA, length(colnames(avg_expr))), colnames(avg_expr))
  }
  roi_vals <- build_roi_values(ct_vals, ct_to_roi_leaf)
  # Get coefficient for title
  coef <- biomarker$mean_coefficient[biomarker$feature == gene]
  direction <- ifelse(coef > 0, "+", "-")
  title <- paste0(gene, " (", direction, abs(round(coef, 2)), ")")
  gene_plots_leaf[[gene]] <- make_heatmap(leaf_map, roi_vals, title,
                                           fill_name = "Expr", fill_scale = "gradient")
}

# Combine into faceted plot using patchwork-like approach with cowplot or gridExtra
# Use patchwork if available, otherwise gridExtra
library(patchwork)

composite1 <- (gene_plots_leaf[[1]] | gene_plots_leaf[[2]]) /
  (gene_plots_leaf[[3]] | gene_plots_leaf[[4]]) +
  plot_annotation(
    title = "LASSO Biomarker Gene Expression on Arabidopsis Leaf Cross-Section",
    subtitle = "Top 4 genes (100% selection frequency) | Coefficient shown in title",
    theme = theme(plot.title = element_text(family = "sans", size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(family = "sans", size = 10, hjust = 0.5, color = "grey30"))
  )

ggsave(paste0(RESULTS, "/figures/ggplantmap_composite1_lasso_leaf.png"),
       plot = composite1, height = 10, width = 12, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_composite1_lasso_leaf.svg"),
       plot = composite1, height = 10, width = 12, bg = "white")
cat("Saved: ggplantmap_composite1_lasso_leaf.png/.svg\n\n")

# ============================================================
# COMPOSITE 2: 4 LASSO genes on root tip cross-section (faceted)
# ============================================================
cat("=== Composite 2: 4 LASSO genes on root tip ===\n")

gene_plots_root <- list()
for (gene in top_genes) {
  if (gene %in% rownames(avg_expr)) {
    ct_vals <- as.numeric(avg_expr[gene, ])
    names(ct_vals) <- colnames(avg_expr)
  } else {
    ct_vals <- setNames(rep(NA, length(colnames(avg_expr))), colnames(avg_expr))
  }
  roi_vals <- build_roi_values(ct_vals, ct_to_roi_root)
  coef <- biomarker$mean_coefficient[biomarker$feature == gene]
  direction <- ifelse(coef > 0, "+", "-")
  title <- paste0(gene, " (", direction, abs(round(coef, 2)), ")")
  gene_plots_root[[gene]] <- make_heatmap(root_map, roi_vals, title,
                                           fill_name = "Expr", fill_scale = "gradient")
}

composite2 <- (gene_plots_root[[1]] | gene_plots_root[[2]]) /
  (gene_plots_root[[3]] | gene_plots_root[[4]]) +
  plot_annotation(
    title = "LASSO Biomarker Gene Expression on Arabidopsis Root Tip Cross-Section",
    subtitle = "Top 4 genes | Cell types mapped: Epidermal->Epidermis, Mesophyll->Cortex, Stele->vascular, Meristematic->Columella",
    theme = theme(plot.title = element_text(family = "sans", size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(family = "sans", size = 9, hjust = 0.5, color = "grey30"))
  )

ggsave(paste0(RESULTS, "/figures/ggplantmap_composite2_lasso_root.png"),
       plot = composite2, height = 10, width = 12, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_composite2_lasso_root.svg"),
       plot = composite2, height = 10, width = 12, bg = "white")
cat("Saved: ggplantmap_composite2_lasso_root.png/.svg\n\n")

# ============================================================
# COMPOSITE 3: CCC communication strength on leaf + root (side-by-side)
# ============================================================
cat("=== Composite 3: CCC communication on leaf + root ===\n")

# Overall communication strength per cell type (row sums = outgoing, col sums = incoming)
# Use total (outgoing + incoming) per cell type
outgoing <- rowSums(ccc)
incoming <- colSums(ccc)
total_comm <- outgoing + incoming
names(total_comm) <- rownames(ccc)

# Leaf
roi_vals_comm_leaf <- build_roi_values(total_comm, ct_to_roi_leaf)
p_comm_leaf <- make_heatmap(leaf_map, roi_vals_comm_leaf,
                             "CCC Communication Strength\n(Leaf)",
                             fill_name = "Strength", fill_scale = "viridis")

# Root
roi_vals_comm_root <- build_roi_values(total_comm, ct_to_roi_root)
p_comm_root <- make_heatmap(root_map, roi_vals_comm_root,
                             "CCC Communication Strength\n(Root Tip)",
                             fill_name = "Strength", fill_scale = "viridis")

composite3 <- (p_comm_leaf | p_comm_root) +
  plot_annotation(
    title = "Cell-Cell Communication Strength on Plant Anatomy",
    subtitle = "Total outgoing + incoming communication per cell type | Ca2+ dominant pathway",
    theme = theme(plot.title = element_text(family = "sans", size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(family = "sans", size = 10, hjust = 0.5, color = "grey30"))
  )

ggsave(paste0(RESULTS, "/figures/ggplantmap_composite3_ccc_communication.png"),
       plot = composite3, height = 6, width = 12, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_composite3_ccc_communication.svg"),
       plot = composite3, height = 6, width = 12, bg = "white")
cat("Saved: ggplantmap_composite3_ccc_communication.png/.svg\n\n")

# ============================================================
# COMPOSITE 4: Ca2+ signaling strength on leaf + root (side-by-side)
# ============================================================
cat("=== Composite 4: Ca2+ signaling on leaf + root ===\n")

# Ca2+ signaling strength per cell type
ca2_signal <- signal_comm %>%
  filter(Signal == "Ca2+") %>%
  group_by(Source) %>%
  summarise(strength = sum(Prob, na.rm = TRUE)) %>%
  pull(strength, name = Source)

# Also add incoming
ca2_incoming <- signal_comm %>%
  filter(Signal == "Ca2+") %>%
  group_by(Target) %>%
  summarise(strength = sum(Prob, na.rm = TRUE)) %>%
  pull(strength, name = Target)

ca2_total <- ca2_signal[rownames(ccc)] + ca2_incoming[rownames(ccc)]
names(ca2_total) <- rownames(ccc)

# Leaf
roi_vals_ca2_leaf <- build_roi_values(ca2_total, ct_to_roi_leaf)
p_ca2_leaf <- make_heatmap(leaf_map, roi_vals_ca2_leaf,
                            "Ca2+ Signaling Strength\n(Leaf)",
                            fill_name = "Ca2+ strength", fill_scale = "viridis")

# Root
roi_vals_ca2_root <- build_roi_values(ca2_total, ct_to_roi_root)
p_ca2_root <- make_heatmap(root_map, roi_vals_ca2_root,
                            "Ca2+ Signaling Strength\n(Root Tip)",
                            fill_name = "Ca2+ strength", fill_scale = "viridis")

composite4 <- (p_ca2_leaf | p_ca2_root) +
  plot_annotation(
    title = "Ca2+ Signaling Strength on Plant Anatomy",
    subtitle = "Dominant CCC pathway (total strength = 249.5, 4.6x next pathway)",
    theme = theme(plot.title = element_text(family = "sans", size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(family = "sans", size = 10, hjust = 0.5, color = "grey30"))
  )

ggsave(paste0(RESULTS, "/figures/ggplantmap_composite4_ca2_signaling.png"),
       plot = composite4, height = 6, width = 12, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_composite4_ca2_signaling.svg"),
       plot = composite4, height = 6, width = 12, bg = "white")
cat("Saved: ggplantmap_composite4_ca2_signaling.png/.svg\n\n")

# ============================================================
# COMPOSITE 5: LASSO panel weighted score on leaf (diverging)
# ============================================================
cat("=== Composite 5: LASSO panel score on leaf ===\n")

# Calculate LASSO panel weighted score per cell type
# score = sum(coef * expr) for all LASSO genes in atlas
lasso_genes_in_atlas <- intersect(rownames(avg_expr), biomarker$feature[biomarker$is_stable])
panel_scores <- rep(0, ncol(avg_expr))
names(panel_scores) <- colnames(avg_expr)
for (g in lasso_genes_in_atlas) {
  coef <- biomarker$mean_coefficient[biomarker$feature == g]
  expr <- as.numeric(avg_expr[g, ])
  panel_scores <- panel_scores + coef * expr
}

cat("LASSO panel scores per cell type:\n")
print(round(panel_scores, 4))

# Leaf (diverging scale - positive=flight, negative=ground)
roi_vals_score_leaf <- build_roi_values(panel_scores, ct_to_roi_leaf)
p_score_leaf <- make_heatmap(leaf_map, roi_vals_score_leaf,
                              "LASSO Panel Weighted Score\n(Leaf)",
                              fill_name = "Score\n(flight+ / ground-)",
                              fill_scale = "diverging")

# Root
roi_vals_score_root <- build_roi_values(panel_scores, ct_to_roi_root)
p_score_root <- make_heatmap(root_map, roi_vals_score_root,
                              "LASSO Panel Weighted Score\n(Root Tip)",
                              fill_name = "Score\n(flight+ / ground-)",
                              fill_scale = "diverging")

composite5 <- (p_score_leaf | p_score_root) +
  plot_annotation(
    title = "LASSO Biomarker Panel Weighted Score on Plant Anatomy",
    subtitle = "Score = sum(coefficient x expression) | Positive (orange) = flight-like, Negative (blue) = ground-like",
    theme = theme(plot.title = element_text(family = "sans", size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(family = "sans", size = 9, hjust = 0.5, color = "grey30"))
  )

ggsave(paste0(RESULTS, "/figures/ggplantmap_composite5_lasso_panel_score.png"),
       plot = composite5, height = 6, width = 12, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ggplantmap_composite5_lasso_panel_score.svg"),
       plot = composite5, height = 6, width = 12, bg = "white")
cat("Saved: ggplantmap_composite5_lasso_panel_score.png/.svg\n\n")

cat("=== All 5 composites generated ===\n")
cat(paste0("Files in ", RESULTS, "/figures/:\n"))
cat("  ggplantmap_composite1_lasso_leaf.png/.svg\n")
cat("  ggplantmap_composite2_lasso_root.png/.svg\n")
cat("  ggplantmap_composite3_ccc_communication.png/.svg\n")
cat("  ggplantmap_composite4_ca2_signaling.png/.svg\n")
cat("  ggplantmap_composite5_lasso_panel_score.png/.svg\n")
