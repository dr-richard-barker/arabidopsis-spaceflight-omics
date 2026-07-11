#!/usr/bin/env Rscript
# ============================================================
# Figure 8 panel a1: CBL9-CIPK23-AKT1 cascade sub-network
# ggraph + igraph + graphlayouts (stress-majorization)
# Nodes: colored by type, sized by guard cell expression
# Edges: Ca2+ (blue), K+ (green), crosstalk CIPK23->AKT1 (orange)
# Key edges labeled with co-expression r values via midpoint annotation
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
library(igraph)
library(ggraph)
library(ggplot2)
library(graphlayouts)
library(dplyr)

# ---- 1. Load data ----
expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"), stringsAsFactors = FALSE)
coexpr <- read.csv(paste0(RESULTS, "/tables/ca2_k_coexpression_Guard.csv"),
                   stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
rownames(coexpr) <- gsub("_Guard", "", rownames(coexpr))

guard_expr <- function(alias) {
  v <- expr$mean_expr[expr$alias == alias & expr$celltype == "Guard"]
  if (length(v) == 0) return(0.01)
  v[1]
}
get_coexpr <- function(g1, g2) {
  if (g1 %in% rownames(coexpr) && g2 %in% colnames(coexpr))
    return(round(coexpr[g1, g2], 3))
  return(NA)
}

# ---- 2. Node definitions ----
nodes <- data.frame(
  name  = c("Ca\u00b2\u207a", "CBL9", "CBL1", "CBL2", "CML24", "CaM3",
            "CIPK23", "CDPK4", "CDPK3", "AKT1", "KC1", "ANNAT1",
            "K\u207a uptake"),
  type  = c("ligand", "sensor_CBL", "sensor_CBL", "sensor_CBL",
            "sensor_CaM", "sensor_CaM",
            "kinase_CIPK", "kinase_CDPK", "kinase_CDPK",
            "channel_K", "channel_K", "target_ca2",
            "output"),
  stringsAsFactors = FALSE
)
nodes$guard_expr_plot <- sapply(nodes$name, function(a) {
  if (grepl("Ca|K\u207a", a)) return(0.12)
  guard_expr(a)
})

cat("Node guard cell expression:\n")
print(nodes[, c("name", "guard_expr_plot")])

# ---- 3. Edge definitions ----
edges <- data.frame(
  from  = c("Ca\u00b2\u207a", "Ca\u00b2\u207a", "Ca\u00b2\u207a", "Ca\u00b2\u207a", "Ca\u00b2\u207a",
            "CBL9", "CBL1", "CBL2",
            "CML24",
            "CaM3", "CaM3",
            "CIPK23",
            "CDPK4",
            "AKT1", "KC1"),
  to    = c("CBL9", "CBL1", "CBL2", "CML24", "CaM3",
            "CIPK23", "CIPK23", "CIPK23",
            "KC1",
            "CDPK4", "CDPK3",
            "AKT1",
            "ANNAT1",
            "K\u207a uptake", "K\u207a uptake"),
  edge_type = c(rep("ca2_binding", 5),
                rep("ca2_activation", 3),
                "ca2_activation",
                rep("ca2_activation", 2),
                "crosstalk",
                "ca2_activation",
                "k_output", "k_output"),
  stringsAsFactors = FALSE
)
edges$coexpr_r <- mapply(get_coexpr, edges$from, edges$to)
edges$edge_label <- ifelse(!is.na(edges$coexpr_r) & abs(edges$coexpr_r) >= 0.1,
                           sprintf("r=%.2f", edges$coexpr_r), "")

cat("\nKey edges with co-expression r:\n")
print(edges[edges$edge_label != "", c("from", "to", "edge_type", "coexpr_r")])

# ---- 4. Build igraph ----
g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

# ---- 5. Stress-majorization layout ----
set.seed(42)
layout_coords <- layout_with_stress(g)
V(g)$x <- layout_coords[, 1]
V(g)$y <- layout_coords[, 2]

# ---- 6. Compute edge midpoints for r-value labels ----
node_xy <- data.frame(
  name = V(g)$name,
  x    = V(g)$x,
  y    = V(g)$y,
  stringsAsFactors = FALSE
)
edge_df <- edges
edge_df$x_from <- node_xy$x[match(edge_df$from, node_xy$name)]
edge_df$y_from <- node_xy$y[match(edge_df$from, node_xy$name)]
edge_df$x_to   <- node_xy$x[match(edge_df$to,   node_xy$name)]
edge_df$y_to   <- node_xy$y[match(edge_df$to,   node_xy$name)]
edge_df$x_mid  <- (edge_df$x_from + edge_df$x_to) / 2
edge_df$y_mid  <- (edge_df$y_from + edge_df$y_to) / 2
label_df <- edge_df[edge_df$edge_label != "", ]

# ---- 7. Color/size scales ----
type_colors <- c(
  "ligand"       = "#E9ED4C",
  "sensor_CBL"   = "#0279EE",
  "sensor_CaM"   = "#FD9BED",
  "kinase_CIPK"  = "#75A025",
  "kinase_CDPK"  = "#FF9400",
  "channel_K"    = "#009E73",
  "target_ca2"   = "#56B4E9",
  "output"       = "#ECE9E2"
)
type_labels <- c(
  "ligand"       = "Ca\u00b2\u207a signal",
  "sensor_CBL"   = "CBL sensor",
  "sensor_CaM"   = "CaM/CML sensor",
  "kinase_CIPK"  = "CIPK kinase",
  "kinase_CDPK"  = "CDPK kinase",
  "channel_K"    = "K\u207a channel",
  "target_ca2"   = "Ca\u00b2\u207a target",
  "output"       = "K\u207a uptake"
)
edge_colors <- c(
  "ca2_binding"    = "#0279EE",
  "ca2_activation" = "#0279EE",
  "crosstalk"      = "#FF9400",
  "k_output"       = "#75A025"
)
edge_widths <- c(
  "ca2_binding"    = 1.0,
  "ca2_activation" = 1.0,
  "crosstalk"      = 2.2,
  "k_output"       = 1.6
)

# ---- 8. Build ggraph ----
p <- ggraph(g, layout = "manual", x = x, y = y) +
  geom_edge_link(
    aes(color = edge_type, width = edge_type),
    arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
    end_cap   = circle(0.42, "cm"),
    start_cap = circle(0.30, "cm"),
    alpha = 0.80
  ) +
  scale_edge_color_manual(values = edge_colors, guide = "none") +
  scale_edge_width_manual(values = edge_widths, guide = "none") +
  # Edge r-value labels at midpoints
  annotate("label",
    x = label_df$x_mid, y = label_df$y_mid,
    label = label_df$edge_label,
    size = 2.8, fontface = "italic", color = "grey20",
    label.padding = unit(0.08, "cm"), label.size = 0,
    fill = alpha("white", 0.80)
  ) +
  # Nodes
  geom_node_point(
    aes(fill = type, size = guard_expr_plot),
    shape = 21, color = "grey20", stroke = 0.8
  ) +
  scale_fill_manual(
    values = type_colors,
    labels = type_labels,
    name   = "Node type",
    guide  = guide_legend(override.aes = list(size = 5))
  ) +
  scale_size_continuous(
    range  = c(4, 13),
    name   = "Guard cell\nexpression",
    breaks = c(0.05, 0.15, 0.35, 0.60),
    labels = c("0.05", "0.15", "0.35", "0.60")
  ) +
  # Node labels
  geom_node_text(
    aes(label = name),
    size = 3.4, fontface = "bold", family = "sans",
    repel = TRUE, max.overlaps = 30,
    bg.colour = "white", bg.r = 0.15,
    box.padding = 0.35, point.padding = 0.45
  ) +
  # Title annotation
  annotate("text",
    x = mean(range(layout_coords[, 1])),
    y = max(layout_coords[, 2]) + 0.65,
    label = "CBL9\u2013CIPK23\u2013AKT1 Cascade",
    size = 5.2, fontface = "bold", family = "sans", hjust = 0.5
  ) +
  annotate("text",
    x = mean(range(layout_coords[, 1])),
    y = max(layout_coords[, 2]) + 0.28,
    label = "Node size \u221d guard cell expression  |  Orange edge = Ca\u00b2\u207a/K\u207a crosstalk",
    size = 3.2, family = "sans", hjust = 0.5, color = "grey30"
  ) +
  theme_void(base_family = "sans") +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(15, 10, 10, 10),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ---- 9. Save ----
out_file <- paste0(WORK, "/figure8a1_cascade.png")
ggsave(out_file, p, width = 9, height = 8, dpi = 300, bg = "white")
cat("\nSaved:", out_file, "-", file.size(out_file), "bytes\n")
