#!/usr/bin/env Rscript
# ============================================================
# ggKEGG Figure 8 Panel A: CBL9-CIPK23-AKT1 cascade on KEGG
# ath04075 (Plant hormone signal transduction)
# with guard cell expression overlay
# ============================================================
# Key design:
#   - KEGG ath04075 PNG as background (annotation_raster)
#   - Highlight KEGG nodes containing circuit genes (CDPK/CIPK, CaM)
#   - Custom overlay nodes for genes NOT in KEGG (CBL9, CBL1, CBL2, CML24, AKT1, KC1, ANNAT1)
#   - Cascade arrows: Ca2+ -> CBL9 -> CIPK23 -> AKT1 -> K+ uptake
#   - Guard cell expression as node fill color (blue gradient)
#   - Co-expression r values on key edges
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

library(ggplot2)
library(png)
library(dplyr)
library(grid)
library(KEGGREST)

# ---- 1. Download KEGG pathway image ----
cat("Downloading KEGG ath04075 image...\n")
# keggGet('ath04075','image') returns a 3D raster array [height, width, 4] directly
kegg_img <- keggGet("ath04075", "image")
tmp_file <- "/tmp/ath04075_map.png"
writePNG(kegg_img, tmp_file)
img_w <- ncol(kegg_img)  # 1095
img_h <- nrow(kegg_img)  # 1726
cat("KEGG image:", img_w, "x", img_h, "\n")

# ---- 2. Guard cell expression data (verified from ca2_k_expr_per_celltype.csv) ----
expr_data <- data.frame(
  gene = c("CBL9", "CBL1", "CBL2", "CML24", "CaM3", "CaM7",
           "CIPK23", "CDPK4", "CDPK3", "AKT1", "KC1", "ANNAT1"),
  guard_expr = c(0.0751, 0.2129, 0.1813, 0.6033, 0.1062, 0.0466,
                 0.1070, 0.0304, 0.0864, 0.1915, 0.5352, 0.0293),
  guard_pct = c(3.4, 11.2, 7.9, 27.0, 5.6, 2.2,
                6.7, 2.2, 5.6, 11.2, 29.2, 2.2),
  in_kegg = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
              TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# KEGG node positions (from KGML, y positive top-to-bottom)
# CDPK/CIPK node: x=465, y=1490 (KO:K13412, "CDPK1...")
# CaM node: x=493, y=1632 (KO:K02183, "CAM4...")
kegg_pos <- data.frame(
  gene = c("CaM3", "CaM7", "CIPK23", "CDPK4", "CDPK3"),
  kegg_x = c(493, 493, 465, 465, 465),
  kegg_y = c(1632, 1632, 1490, 1490, 1490),
  stringsAsFactors = FALSE
)

# Custom positions for non-KEGG genes — arranged as cascade
# CBLs: upstream (above CDPK/CIPK node, left side)
# CML24: near CaM node (right side)
# AKT1, KC1, ANNAT1: downstream (below CDPK/CIPK node)
custom_pos <- data.frame(
  gene = c("CBL9", "CBL1", "CBL2", "CML24", "AKT1", "KC1", "ANNAT1"),
  x = c(180, 280, 380, 680, 180, 330, 480),
  y = c(1380, 1330, 1280, 1632, 1720, 1770, 1820),
  stringsAsFactors = FALSE
)

expr_data <- expr_data %>%
  left_join(kegg_pos, by = "gene") %>%
  left_join(custom_pos, by = "gene") %>%
  mutate(
    plot_x = ifelse(in_kegg, kegg_x, x),
    plot_y = ifelse(in_kegg, kegg_y, y)
  )

# ---- 3. Co-expression r values (verified from ca2_k_coexpression_Guard.csv) ----
coexp_edges <- data.frame(
  from = c("CBL9", "CBL9", "CIPK23", "CDPK4", "KC1", "ANNAT1"),
  to = c("CIPK23", "AKT1", "CaM3", "ANNAT1", "CaM3", "CaM3"),
  r = c(-0.047, 0.337, 0.460, 0.559, 0.237, 0.227),
  stringsAsFactors = FALSE
)

# Get positions for edge labels
get_pos <- function(gene) {
  row <- expr_data[expr_data$gene == gene, ]
  c(row$plot_x, row$plot_y)
}

edge_labels <- data.frame()
for (i in seq_len(nrow(coexp_edges))) {
  from_pos <- get_pos(coexp_edges$from[i])
  to_pos <- get_pos(coexp_edges$to[i])
  edge_labels <- rbind(edge_labels, data.frame(
    x = (from_pos[1] + to_pos[1]) / 2,
    y = (from_pos[2] + to_pos[2]) / 2,
    label = sprintf("r=%.2f", coexp_edges$r[i])
  ))
}

# ---- 4. Cascade arrows ----
# Main cascade: Ca2+ -> CBL9 -> CIPK23 -> AKT1 -> K+ uptake
# Also: CBL9 -> AKT1 (direct co-expression edge)
cascade <- data.frame(
  x_from = c(300, 180, 465, 180, 180),
  y_from = c(1200, 1380, 1490, 1380, 1720),
  x_to = c(180, 465, 180, 180, 100),
  y_to = c(1380, 1490, 1720, 1720, 1900),
  color = c("#E9ED4C", "#0279EE", "#75A025", "#0279EE", "#FF9400"),
  linewidth = c(2.0, 2.0, 2.5, 1.5, 2.5),
  label = c("", "r=-0.05", "", "r=0.34", ""),
  stringsAsFactors = FALSE
)

# ---- 5. Build the plot ----
cat("Building ggplot...\n")

# Non-KEGG genes for custom nodes
non_kegg <- expr_data[!expr_data$in_kegg, ]

p <- ggplot() +
  # KEGG background image
  annotation_raster(kegg_img, xmin = 0, xmax = img_w,
                    ymin = 0, ymax = img_h, interpolate = TRUE) +

  # --- Cascade arrows (drawn before nodes so nodes sit on top) ---
  geom_segment(data = cascade,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = cascade$color, linewidth = cascade$linewidth,
               alpha = 0.85, show.legend = FALSE) +

  # Co-expression edge labels on cascade arrows
  geom_text(data = edge_labels,
            aes(x = x, y = y, label = label),
            size = 3, fontface = "italic", color = "grey20",
            nudge_x = 20, nudge_y = 15) +

  # --- Ca2+ source node (yellow) ---
  annotate("point", x = 300, y = 1200, fill = "#E9ED4C",
           shape = 21, size = 12, stroke = 2, color = "grey30") +
  annotate("text", x = 300, y = 1200 - 35, label = "Ca\u00b2\u207a",
           size = 5, fontface = "bold") +

  # --- K+ uptake target node (green) ---
  annotate("point", x = 100, y = 1900, fill = "#75A025",
           shape = 21, size = 12, stroke = 2, color = "grey30") +
  annotate("text", x = 100, y = 1900 + 35, label = "K\u207a uptake",
           size = 4.5, fontface = "bold") +

  # --- Custom nodes for non-KEGG genes (filled by guard cell expression) ---
  geom_point(data = non_kegg,
             aes(x = plot_x, y = plot_y, fill = guard_expr),
             shape = 21, size = 8, stroke = 1.5, color = "grey30") +
  scale_fill_gradient(low = "white", high = "#0279EE",
                      name = "Guard cell\nexpression",
                      limits = c(0, 0.65),
                      breaks = c(0, 0.2, 0.4, 0.6),
                      guide = guide_colorbar(
                        title.vjust = 0.8,
                        barwidth = unit(0.5, "cm"),
                        barheight = unit(3, "cm")
                      )) +

  # Labels for custom nodes
  geom_text(data = non_kegg,
            aes(x = plot_x, y = plot_y - 25,
                label = paste0(gene, "\n(", round(guard_expr, 3), ")")),
            size = 2.8, fontface = "bold", lineheight = 0.85) +

  # --- Highlight CDPK/CIPK KEGG node (orange) ---
  annotate("rect", xmin = 465 - 75, xmax = 465 + 75,
           ymin = 1490 - 20, ymax = 1490 + 20,
           fill = "#FF9400", alpha = 0.35, color = "#FF6600", linewidth = 2) +
  annotate("text", x = 465, y = 1490 - 32,
           label = "CDPK4/3/CIPK23\n(guard: 0.075)",
           size = 2.8, fontface = "bold", color = "#CC4400", lineheight = 0.85) +

  # --- Highlight CaM KEGG node (pink) ---
  annotate("rect", xmin = 493 - 75, xmax = 493 + 75,
           ymin = 1632 - 20, ymax = 1632 + 20,
           fill = "#FD9BED", alpha = 0.35, color = "#DD7BCD", linewidth = 2) +
  annotate("text", x = 493, y = 1632 - 32,
           label = "CaM3/7\n(guard: 0.076)",
           size = 2.8, fontface = "bold", color = "#AA3399", lineheight = 0.85) +

  # --- Title and subtitle ---
  annotate("text", x = 547, y = 60,
           label = "CBL9-CIPK23-AKT1 Cascade on KEGG Plant Hormone Signaling",
           size = 6, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 547, y = 95,
           label = "Guard cell expression overlay  |  ath04075",
           size = 4, hjust = 0.5, color = "grey40") +

  # --- Coordinate system ---
  coord_cartesian(xlim = c(0, img_w), ylim = c(0, img_h), expand = FALSE) +
  scale_y_reverse() +

  # --- Theme ---
  theme_void() +
  theme(
    legend.position = c(0.92, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- 6. Save ----
cat("Saving figure...\n")
out_file <- paste0(WORK, "/figure8a_ggkegg.png")
ggsave(out_file, p, width = 10, height = 15, dpi = 300, bg = "white")
cat("Saved:", out_file, "-", file.size(out_file), "bytes\n")

# Also save SVG version
out_svg <- paste0(WORK, "/figure8a_ggkegg.svg")
# SVG won't embed the raster well, so we skip SVG for the KEGG panel
# (the PNG is the primary deliverable)
cat("Done.\n")
