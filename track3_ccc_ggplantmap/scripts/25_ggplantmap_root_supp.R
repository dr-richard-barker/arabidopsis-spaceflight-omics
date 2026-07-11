#!/usr/bin/env Rscript
# ============================================================
# Supplementary Figure: Root mature cross-section
# ggPm.At.rootmatur.crosssection — 8 ROIs
# ROI -> cell type mapping:
#   Atrichoblast  -> Epidermal
#   Trichoblast   -> Epidermal
#   Cortex        -> Mesophyll
#   Endodermis    -> Mesophyll
#   Pericycle     -> Stele
#   Phloem        -> Stele
#   Xylem         -> Stele
#   Procambium    -> Meristematic
# Three sub-maps: Ca2+ CCC strength, AKT1 expression, CBL9 expression
# ============================================================
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
library(patchwork)

# ---- 1. Inspect root map ----
root_map <- ggPm.At.rootmatur.crosssection
cat("Root mature cross-section ROIs:\n")
print(unique(root_map[, c("ROI.name", "ROI.id")]))

# ---- 2. ROI -> cell type mapping ----
roi_to_celltype <- c(
  "Atrichoblast" = "Epidermal",
  "Trichoblast"  = "Epidermal",
  "Cortex"       = "Mesophyll",
  "Endodermis"   = "Mesophyll",
  "Pericycle"    = "Stele",
  "Phloem"       = "Stele",
  "Xylem"        = "Stele",
  "Procambium"   = "Meristematic"
)

# ---- 3. Data values ----
# Ca2+ CCC outgoing strength
ca2_strength <- c(Guard=57.12, Meristematic=50.98, Epidermal=49.57, Mesophyll=47.49, Stele=44.31)
# K+ CCC outgoing strength
k_strength   <- c(Guard=4.35, Meristematic=4.49, Epidermal=6.06, Mesophyll=3.76, Stele=5.84)

# Gene expression
expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"), stringsAsFactors = FALSE)
get_expr <- function(alias, celltype) {
  v <- expr$mean_expr[expr$alias == alias & expr$celltype == celltype]
  if (length(v) == 0) return(NA)
  v[1]
}

# ---- 4. Build value tables ----
roi_names <- unique(root_map$ROI.name)

make_vals <- function(value_vec) {
  data.frame(
    ROI.name = roi_names,
    value    = value_vec[roi_to_celltype[roi_names]],
    stringsAsFactors = FALSE
  )
}

make_gene_vals <- function(gene_alias) {
  data.frame(
    ROI.name = roi_names,
    value    = sapply(roi_to_celltype[roi_names], function(ct) get_expr(gene_alias, ct)),
    stringsAsFactors = FALSE
  )
}

ca2_vals  <- make_vals(ca2_strength)
k_vals    <- make_vals(k_strength)
akt1_vals <- make_gene_vals("AKT1")
cbl9_vals <- make_gene_vals("CBL9")

cat("\nCa2+ CCC strength by ROI:\n"); print(ca2_vals)
cat("\nAKT1 expression by ROI:\n");   print(akt1_vals)
cat("\nCBL9 expression by ROI:\n");   print(cbl9_vals)

# ---- 5. Merge ----
root_ca2  <- ggPlantmap.merge(root_map, ca2_vals,  id.x = "ROI.name", id.y = "ROI.name")
root_k    <- ggPlantmap.merge(root_map, k_vals,    id.x = "ROI.name", id.y = "ROI.name")
root_akt1 <- ggPlantmap.merge(root_map, akt1_vals, id.x = "ROI.name", id.y = "ROI.name")
root_cbl9 <- ggPlantmap.merge(root_map, cbl9_vals, id.x = "ROI.name", id.y = "ROI.name")

# ---- 6. Plot function ----
plot_root <- function(merged_data, color_high, legend_name, title_label, subtitle_label,
                      vmin = NULL, vmax = NULL) {
  p <- ggPlantmap.heatmap(
    map.quant   = merged_data,
    value.quant = value,
    show.legend = TRUE
  ) +
    scale_fill_gradient(
      low    = "#ECE9E2",
      high   = color_high,
      name   = legend_name,
      limits = if (!is.null(vmin)) c(vmin, vmax) else NULL,
      na.value = "grey80"
    ) +
    labs(title = title_label, subtitle = subtitle_label) +
    theme_void(base_family = "sans") +
    theme(
      plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
      legend.title  = element_text(size = 8, face = "bold"),
      legend.text   = element_text(size = 7),
      plot.margin   = margin(3, 3, 3, 3)
    )
  p
}

# ---- 7. Build 4 sub-maps ----
p_ca2  <- plot_root(root_ca2,  "#0279EE",
                    "Ca\u00b2\u207a CCC\nstrength",
                    "(i) Ca\u00b2\u207a CCC outgoing strength",
                    "Meristematic highest (51.0)",
                    vmin = 40, vmax = 55)

p_k    <- plot_root(root_k,    "#75A025",
                    "K\u207a CCC\nstrength",
                    "(ii) K\u207a CCC outgoing strength",
                    "Epidermal highest (6.1)",
                    vmin = 3.5, vmax = 6.5)

p_akt1 <- plot_root(root_akt1, "#FF9400",
                    "AKT1\nexpr",
                    "(iii) AKT1 expression",
                    "Epidermal highest (0.367)",
                    vmin = 0.05, vmax = 0.40)

p_cbl9 <- plot_root(root_cbl9, "#FD9BED",
                    "CBL9\nexpr",
                    "(iv) CBL9 expression",
                    "Meristematic highest",
                    vmin = 0.01, vmax = 0.15)

# ---- 8. Combine ----
p_combined <- p_ca2 + p_k + p_akt1 + p_cbl9 +
  plot_layout(ncol = 4) +
  plot_annotation(
    title    = "Root mature cross-section: Ca\u00b2\u207a/K\u207a circuit",
    subtitle = "CCC signaling strength and cascade gene expression by root cell type",
    caption  = paste0("ROI mapping: Atrichoblast/Trichoblast \u2192 Epidermal; ",
                      "Cortex/Endodermis \u2192 Mesophyll; ",
                      "Pericycle/Phloem/Xylem \u2192 Stele; ",
                      "Procambium \u2192 Meristematic"),
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5, family = "sans"),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40", family = "sans"),
      plot.caption  = element_text(size = 7,  hjust = 0.5, color = "grey55", family = "sans")
    )
  )

# ---- 9. Save ----
out_file <- paste0(WORK, "/SuppFig_root_crosssection.png")
ggsave(out_file, p_combined, width = 16, height = 5.5, dpi = 300, bg = "white")
cat("\nSaved:", out_file, "-", file.size(out_file), "bytes\n")
