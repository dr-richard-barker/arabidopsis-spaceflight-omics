#!/usr/bin/env Rscript
# ============================================================
# Supplementary Figure: Seed developmental series
# ggPm.At.seed.devseries — 5 stages x 3 genes (CBL9, CIPK23, AKT1)
# Proxy cell-type mapping:
#   Embryo Proper        -> Meristematic (+ Guard avg)
#   Micropylar Endosperm -> Stele
#   Peripheral Endosperm -> Stele
#   Chalazal Endosperm   -> Mesophyll
#   Distal Seed Coat     -> Epidermal
#   Chalazal Seed Coat   -> Epidermal
# Stages: Preglobular, Globular, Heart, Linearcotyledon, Maturegreen
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

# ---- 1. Inspect seed map ----
seed_map <- ggPm.At.seed.devseries
cat("Seed map ROI names (unique):\n")
print(unique(seed_map$ROI.name))
cat("\nStages:\n")
print(unique(seed_map$Stage))
cat("\nParts:\n")
print(unique(seed_map$Part))
cat("\nRegions:\n")
print(unique(seed_map$Region))

# ---- 2. Expression data ----
expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"), stringsAsFactors = FALSE)

get_expr <- function(alias, celltype) {
  v <- expr$mean_expr[expr$alias == alias & expr$celltype == celltype]
  if (length(v) == 0) return(NA)
  v[1]
}

# Proxy mapping: seed ROI Region -> cell type
# Based on biological correspondence
region_to_celltype <- c(
  "Embryo Proper"          = "Meristematic",
  "Micropylar Endosperm"   = "Stele",
  "Peripheral Endosperm"   = "Stele",
  "Chalazal Endosperm"     = "Mesophyll",
  "Distal Seed Coat"       = "Epidermal",
  "Chalazal Seed Coat"     = "Epidermal"
)

# Gene expression by cell type
genes <- c("CBL9", "CIPK23", "AKT1")
gene_expr <- lapply(genes, function(g) {
  sapply(c("Meristematic", "Stele", "Mesophyll", "Epidermal", "Guard"), function(ct) {
    get_expr(g, ct)
  })
})
names(gene_expr) <- genes
cat("\nGene expression by cell type:\n")
print(do.call(rbind, gene_expr))

# ---- 3. Build value tables for each gene ----
# Get unique ROI names from seed map
roi_names <- unique(seed_map$ROI.name)
cat("\nAll ROI names:\n")
print(roi_names)

# Parse ROI name: format is "Stage.part.Region" e.g. "Globular.embryo.Embryo Proper"
# Extract Region from ROI name
parse_region <- function(roi_name) {
  # ROI.name format: Stage.Part.Region (Region may have spaces)
  parts <- strsplit(roi_name, "\\.")[[1]]
  if (length(parts) >= 3) {
    return(paste(parts[3:length(parts)], collapse = "."))
  }
  return(NA)
}

# Actually check the Region column directly
roi_region_map <- seed_map %>%
  distinct(ROI.name, Stage, Part, Region) %>%
  mutate(celltype = region_to_celltype[Region])

cat("\nROI -> celltype mapping:\n")
print(roi_region_map)

# ---- 4. Build merged data for each gene ----
make_gene_vals <- function(gene_alias) {
  roi_region_map %>%
    mutate(value = sapply(celltype, function(ct) {
      if (is.na(ct)) return(NA)
      get_expr(gene_alias, ct)
    })) %>%
    select(ROI.name, value)
}

# ---- 5. Plot function ----
plot_gene_seed <- function(gene_alias, color_high, title_label) {
  vals <- make_gene_vals(gene_alias)
  merged <- ggPlantmap.merge(seed_map, vals, id.x = "ROI.name", id.y = "ROI.name")
  
  # Get expression range for this gene across all cell types
  all_vals <- vals$value[!is.na(vals$value)]
  vmin <- min(all_vals) * 0.9
  vmax <- max(all_vals) * 1.05
  
  ggPlantmap.heatmap(
    map.quant   = merged,
    value.quant = value,
    show.legend = TRUE
  ) +
    scale_fill_gradient(
      low    = "#ECE9E2",
      high   = color_high,
      name   = paste0(gene_alias, "\nexpr"),
      limits = c(vmin, vmax),
      na.value = "grey80"
    ) +
    labs(title = title_label) +
    theme_void(base_family = "sans") +
    theme(
      plot.title   = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text  = element_text(size = 7),
      plot.margin  = margin(3, 3, 3, 3)
    )
}

# ---- 6. Build 3-gene plots ----
p_cbl9   <- plot_gene_seed("CBL9",   "#0279EE", "CBL9 (Ca\u00b2\u207a sensor)")
p_cipk23 <- plot_gene_seed("CIPK23", "#75A025", "CIPK23 (kinase)")
p_akt1   <- plot_gene_seed("AKT1",   "#FF9400", "AKT1 (K\u207a channel)")

# ---- 7. Combine ----
p_combined <- p_cbl9 + p_cipk23 + p_akt1 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title    = "Seed developmental series: CBL9\u2013CIPK23\u2013AKT1 cascade expression",
    subtitle = "5 stages (Preglobular \u2192 Maturegreen) \u00b7 proxy cell-type mapping from scRNA-seq",
    caption  = "Regions: Embryo Proper \u2192 Meristematic; Endosperm \u2192 Stele/Mesophyll; Seed Coat \u2192 Epidermal",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5, family = "sans"),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40", family = "sans"),
      plot.caption  = element_text(size = 7,  hjust = 0.5, color = "grey55", family = "sans")
    )
  )

# ---- 8. Save ----
out_file <- paste0(WORK, "/SuppFig_seed_devseries.png")
ggsave(out_file, p_combined, width = 14, height = 6, dpi = 300, bg = "white")
cat("\nSaved:", out_file, "-", file.size(out_file), "bytes\n")
