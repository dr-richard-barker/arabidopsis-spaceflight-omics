#!/usr/bin/env Rscript
# ============================================================
# Figure 8 panel a2: ggPlantmap leaf cross-section
# Two sub-maps side by side:
#   (i)  Ca2+ CCC outgoing strength  — guard cell dominance
#   (ii) AKT1 expression             — K+ uptake terminus
# Uses ggPlantmap.merge + ggPlantmap.heatmap (continuous fill)
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

# ---- 1. Inspect leaf cross-section ROIs ----
leaf_map <- ggPm.At.leaf.crosssection
cat("Leaf cross-section ROI names:\n")
print(unique(leaf_map$ROI.name))

# ---- 2. ROI -> cell type mapping ----
roi_to_celltype <- c(
  "epidermis.stomata"          = "Guard",
  "Parenchima.palisade"        = "Mesophyll",
  "Parenchima.sponge"          = "Mesophyll",
  "epidermis.adaxial"          = "Epidermal",
  "epidermis.abaxial"          = "Epidermal",
  "vascularbundle.xylem"       = "Stele",
  "vascularbundle.phloem"      = "Stele",
  "vascularbundle.bundlesheet" = "Stele"
)

# ---- 3. Data values ----
ca2_strength <- c(Guard=57.12, Mesophyll=47.49, Epidermal=49.57, Stele=44.31)
akt1_expr    <- c(Guard=0.1915, Mesophyll=0.1970, Epidermal=0.3670, Stele=0.2880)

# ---- 4. Build value tables for ggPlantmap.merge ----
# Need a data.frame with ROI.name and value columns
roi_names <- unique(leaf_map$ROI.name)

ca2_vals <- data.frame(
  ROI.name = roi_names,
  ca2_strength = ca2_strength[roi_to_celltype[roi_names]],
  stringsAsFactors = FALSE
)
akt1_vals <- data.frame(
  ROI.name = roi_names,
  akt1_expr = akt1_expr[roi_to_celltype[roi_names]],
  stringsAsFactors = FALSE
)

cat("\nCa2+ CCC strength by ROI:\n")
print(ca2_vals)
cat("\nAKT1 expression by ROI:\n")
print(akt1_vals)

# ---- 5. Merge with map ----
leaf_ca2  <- ggPlantmap.merge(leaf_map, ca2_vals,  id.x = "ROI.name", id.y = "ROI.name")
leaf_akt1 <- ggPlantmap.merge(leaf_map, akt1_vals, id.x = "ROI.name", id.y = "ROI.name")

cat("\nMerged Ca2+ columns:", names(leaf_ca2), "\n")

# ---- 6. Sub-map i: Ca2+ CCC outgoing strength ----
p_ca2 <- ggPlantmap.heatmap(
  map.quant   = leaf_ca2,
  value.quant = ca2_strength,
  show.legend = TRUE
) +
  scale_fill_gradient(
    low    = "#ECE9E2",
    high   = "#0279EE",
    name   = "Ca\u00b2\u207a CCC\noutgoing\nstrength",
    limits = c(40, 60),
    breaks = c(44, 48, 52, 57),
    labels = c("44", "48", "52", "57")
  ) +
  labs(
    title    = "(i) Ca\u00b2\u207a CCC outgoing strength",
    subtitle = "Guard cells dominant (57.1)"
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "#0279EE"),
    legend.title  = element_text(size = 8,  face = "bold"),
    legend.text   = element_text(size = 7),
    plot.margin   = margin(5, 5, 5, 5)
  )

# ---- 7. Sub-map ii: AKT1 expression ----
p_akt1 <- ggPlantmap.heatmap(
  map.quant   = leaf_akt1,
  value.quant = akt1_expr,
  show.legend = TRUE
) +
  scale_fill_gradient(
    low    = "#ECE9E2",
    high   = "#75A025",
    name   = "AKT1\nexpression\n(mean)",
    limits = c(0.15, 0.40),
    breaks = c(0.19, 0.25, 0.30, 0.37),
    labels = c("0.19", "0.25", "0.30", "0.37")
  ) +
  labs(
    title    = "(ii) AKT1 expression",
    subtitle = "K\u207a uptake terminus (Epidermal highest)"
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "#75A025"),
    legend.title  = element_text(size = 8,  face = "bold"),
    legend.text   = element_text(size = 7),
    plot.margin   = margin(5, 5, 5, 5)
  )

# ---- 8. Combine with patchwork ----
p_combined <- p_ca2 + p_akt1 +
  plot_annotation(
    title   = "Leaf cross-section: Ca\u00b2\u207a/K\u207a circuit",
    caption = "ROIs: Guard (stomata), Mesophyll (palisade+sponge), Epidermal (adaxial+abaxial), Stele (vascular bundle)",
    theme   = theme(
      plot.title   = element_text(size = 13, face = "bold", hjust = 0.5, family = "sans"),
      plot.caption = element_text(size = 7,  hjust = 0.5, color = "grey50", family = "sans")
    )
  )

# ---- 9. Save ----
out_file <- paste0(WORK, "/figure8a2_leaf.png")
ggsave(out_file, p_combined, width = 10, height = 5.5, dpi = 300, bg = "white")
cat("\nSaved:", out_file, "-", file.size(out_file), "bytes\n")
