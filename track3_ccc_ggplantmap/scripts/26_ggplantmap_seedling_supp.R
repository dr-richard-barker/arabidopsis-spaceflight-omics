#!/usr/bin/env Rscript
# ============================================================
# Supplementary Figure: Young seedling (redesigned)
# ggPm.At.seedling.saltdrought — all 12 ROIs (4 organs x 3 treatments)
# Organ -> cell type proxy:
#   Cotyledon -> avg(Mesophyll, Guard)
#   Hypocotyl -> avg(Epidermal, Stele)
#   Root      -> avg(Stele, Meristematic)
#   Hook      -> avg(all 5 cell types)
# Four sub-maps: Ca2+ CCC strength, K+ CCC strength, AKT1, CBL9
# Note: seedling map has no 'point' column; use direct left_join
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

# ---- 1. Load seedling map ----
seedling_map <- ggPm.At.seedling.saltdrought
cat("Seedling map columns:", names(seedling_map), "\n")
cat("ROI names:\n")
print(unique(seedling_map$ROI.name))

# ---- 2. Data values ----
ca2_strength <- c(Guard=57.12, Meristematic=50.98, Epidermal=49.57, Mesophyll=47.49, Stele=44.31)
k_strength   <- c(Guard=4.35,  Meristematic=4.49,  Epidermal=6.06,  Mesophyll=3.76,  Stele=5.84)

expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"), stringsAsFactors = FALSE)
get_expr <- function(alias, celltype) {
  v <- expr$mean_expr[expr$alias == alias & expr$celltype == celltype]
  if (length(v) == 0) return(NA)
  v[1]
}

# Organ -> value (averaged over proxy cell types)
organ_ca2 <- c(
  Cotyledon = mean(c(ca2_strength["Mesophyll"], ca2_strength["Guard"])),
  Hypocotyl = mean(c(ca2_strength["Epidermal"], ca2_strength["Stele"])),
  Root      = mean(c(ca2_strength["Stele"],     ca2_strength["Meristematic"])),
  Hook      = mean(ca2_strength)
)
organ_k <- c(
  Cotyledon = mean(c(k_strength["Mesophyll"], k_strength["Guard"])),
  Hypocotyl = mean(c(k_strength["Epidermal"], k_strength["Stele"])),
  Root      = mean(c(k_strength["Stele"],     k_strength["Meristematic"])),
  Hook      = mean(k_strength)
)
organ_akt1 <- c(
  Cotyledon = mean(c(get_expr("AKT1","Mesophyll"), get_expr("AKT1","Guard"))),
  Hypocotyl = mean(c(get_expr("AKT1","Epidermal"), get_expr("AKT1","Stele"))),
  Root      = mean(c(get_expr("AKT1","Stele"),     get_expr("AKT1","Meristematic"))),
  Hook      = mean(c(get_expr("AKT1","Mesophyll"), get_expr("AKT1","Guard"),
                     get_expr("AKT1","Epidermal"), get_expr("AKT1","Stele"),
                     get_expr("AKT1","Meristematic")))
)
organ_cbl9 <- c(
  Cotyledon = mean(c(get_expr("CBL9","Mesophyll"), get_expr("CBL9","Guard"))),
  Hypocotyl = mean(c(get_expr("CBL9","Epidermal"), get_expr("CBL9","Stele"))),
  Root      = mean(c(get_expr("CBL9","Stele"),     get_expr("CBL9","Meristematic"))),
  Hook      = mean(c(get_expr("CBL9","Mesophyll"), get_expr("CBL9","Guard"),
                     get_expr("CBL9","Epidermal"), get_expr("CBL9","Stele"),
                     get_expr("CBL9","Meristematic")))
)

cat("\nOrgan Ca2+ CCC strength:\n"); print(round(organ_ca2, 2))
cat("\nOrgan K+ CCC strength:\n");   print(round(organ_k, 2))
cat("\nOrgan AKT1 expression:\n");   print(round(organ_akt1, 4))
cat("\nOrgan CBL9 expression:\n");   print(round(organ_cbl9, 4))

# ---- 3. Build value table (ROI.name -> value via organ column) ----
# The seedling map has an 'organ' column directly
make_seedling_merged <- function(organ_vals) {
  organ_df <- data.frame(
    organ = names(organ_vals),
    value = as.numeric(organ_vals),
    stringsAsFactors = FALSE
  )
  # Join by organ column
  merged <- seedling_map %>%
    left_join(organ_df, by = "organ")
  merged
}

# ---- 4. Plot function using ggPlantmap.heatmap directly ----
# (seedling map has no 'point' column so ggPlantmap.merge won't work)
plot_seedling <- function(organ_vals, color_high, legend_name, title_label, subtitle_label,
                          vmin = NULL, vmax = NULL) {
  merged <- make_seedling_merged(organ_vals)
  
  ar <- (max(merged$y) - min(merged$y)) / (max(merged$x) - min(merged$x))
  
  ggplot(merged, aes(x = x, y = y)) +
    geom_polygon(aes(group = ROI.id, fill = value),
                 colour = "black", linewidth = 0.5, show.legend = TRUE) +
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
      aspect.ratio  = ar,
      plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey40"),
      legend.title  = element_text(size = 8, face = "bold"),
      legend.text   = element_text(size = 7),
      plot.margin   = margin(3, 3, 3, 3)
    )
}

# ---- 5. Build 4 sub-maps ----
p_ca2  <- plot_seedling(organ_ca2,  "#0279EE",
                        "Ca\u00b2\u207a CCC\nstrength",
                        "(i) Ca\u00b2\u207a CCC outgoing strength",
                        "Cotyledon highest (guard+mesophyll avg)",
                        vmin = 44, vmax = 55)

p_k    <- plot_seedling(organ_k,    "#75A025",
                        "K\u207a CCC\nstrength",
                        "(ii) K\u207a CCC outgoing strength",
                        "Hypocotyl highest (epidermal+stele avg)",
                        vmin = 3.5, vmax = 6.5)

p_akt1 <- plot_seedling(organ_akt1, "#FF9400",
                        "AKT1\nexpr",
                        "(iii) AKT1 expression",
                        "Hypocotyl highest (epidermal+stele avg)",
                        vmin = 0.15, vmax = 0.35)

p_cbl9 <- plot_seedling(organ_cbl9, "#FD9BED",
                        "CBL9\nexpr",
                        "(iv) CBL9 expression",
                        "Root highest (stele+meristematic avg)",
                        vmin = 0.05, vmax = 0.15)

# ---- 6. Combine ----
p_combined <- p_ca2 + p_k + p_akt1 + p_cbl9 +
  plot_layout(ncol = 4) +
  plot_annotation(
    title    = "Young seedling: Ca\u00b2\u207a/K\u207a circuit (all treatments shown)",
    subtitle = "CCC signaling strength and cascade gene expression by organ (proxy cell-type mapping)",
    caption  = paste0("Organ proxy: Cotyledon \u2192 Mesophyll+Guard avg; ",
                      "Hypocotyl \u2192 Epidermal+Stele avg; ",
                      "Root \u2192 Stele+Meristematic avg; ",
                      "Hook \u2192 all cell types avg\n",
                      "Treatments shown: control (left), NaCl (center), Sorbitol (right) within each organ"),
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5, family = "sans"),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40", family = "sans"),
      plot.caption  = element_text(size = 7,  hjust = 0.5, color = "grey55", family = "sans")
    )
  )

# ---- 7. Save ----
out_file <- paste0(WORK, "/SuppFig_seedling.png")
ggsave(out_file, p_combined, width = 16, height = 6, dpi = 300, bg = "white")
cat("\nSaved:", out_file, "-", file.size(out_file), "bytes\n")
