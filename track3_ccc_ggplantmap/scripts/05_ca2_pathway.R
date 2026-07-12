#!/usr/bin/env Rscript
# Ca2+ signaling pathway — REFINED network diagram
# Expanded to ~20 nodes with stress-majorization layout (graphlayouts package)
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
suppressMessages({library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(tibble);library(purrr)})
library(igraph)
library(ggraph)
library(graphlayouts)
library(viridis)

cat("=== Ca2+ Pathway Diagram (Refined, Stress-Majorization Layout) ===\n\n")

# Load data
ca2_lr <- read.csv(paste0(RESULTS, "/tables/ccc_lr_pairs.csv"), stringsAsFactors = FALSE) %>%
  filter(Signal == "Ca2+")
ca2_full <- read.csv(paste0(RESULTS, "/tables/ccc_full_communication_table.csv"), stringsAsFactors = FALSE) %>%
  filter(Signal == "Ca2+")
ca2_strength <- ca2_full %>%
  group_by(Interaction_name) %>%
  summarise(total_prob = sum(Prob, na.rm = TRUE)) %>%
  arrange(desc(total_prob))
ca2_expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

cat("Ca2+ LR pairs in DB:", nrow(ca2_lr), "\n")
cat("Measured Ca2+ interactions:", nrow(ca2_strength), "\n\n")

# ============================================================
# Node definitions (20 nodes total)
# ============================================================
node_defs <- data.frame(
  name = c("Ca2+",
           # Ca2+ sensors (CBLs)
           "AT4G17615", "AT5G47100", "AT5G55990",
           # Ca2+ sensor kinases (CIPKs)
           "AT1G30270", "AT3G17510", "AT1G01140",
           # CaM/CML sensors
           "AT1G76650", "AT4G20780", "AT5G37770", "AT3G43810", "AT3G56800",
           # CDPKs
           "AT4G23650", "AT4G09570",
           # Downstream targets
           "AT3G51860", "AT1G08090", "AT2G01980", "AT1G35720",
           # LASSO genes
           "AT4G10250", "AT4G01390"),
  alias = c("Ca2+",
            "CBL1", "CBL9", "CBL2",
            "CIPK23", "CIPK1", "CIPK9",
            "CML38", "CML42", "CML24", "CaM7", "CaM3",
            "CDPK3", "CDPK4",
            "CAX3", "NRT2.1", "SOS1", "ANNAT1",
            "HSP22.0\n(LASSO+)", "TRAF-like\n(LASSO-)"),
  type = c("ligand",
           rep("sensor_CBL", 3),
           rep("kinase_CIPK", 3),
           rep("sensor_CaM", 5),
           rep("kinase_CDPK", 2),
           rep("target", 4),
           rep("LASSO_gene", 2)),
  stringsAsFactors = FALSE
)

cat("Nodes:", nrow(node_defs), "\n")

# ============================================================
# Edge table
# ============================================================
# 1. Primary edges: Ca2+ -> receptors (from measured CCC data)
# Get top interactions involving our node set
node_genes <- setdiff(node_defs$name, "Ca2+")
primary_edges <- data.frame()
for (i in 1:nrow(ca2_lr)) {
  rec <- ca2_lr$Receptor[i]
  lig <- ca2_lr$Ligand[i]
  interaction_name <- ca2_lr$Interaction_name[i]
  strength_val <- ca2_strength$total_prob[ca2_strength$Interaction_name == interaction_name]
  if (length(strength_val) == 0) next

  # Ca2+ -> receptor (if receptor is in our node set)
  if (rec %in% node_genes) {
    primary_edges <- rbind(primary_edges, data.frame(
      from = "Ca2+", to = rec, label = "binds",
      strength = strength_val, interaction = interaction_name,
      edge_type = "primary", stringsAsFactors = FALSE
    ))
  }
  # Ligand -> receptor (if both in node set, ligand is a Ca2+ sensor)
  if (lig %in% node_genes && rec %in% node_genes && lig != rec) {
    primary_edges <- rbind(primary_edges, data.frame(
      from = lig, to = rec, label = "interacts",
      strength = strength_val, interaction = interaction_name,
      edge_type = "primary", stringsAsFactors = FALSE
    ))
  }
}
# Deduplicate and keep top edges
primary_edges <- primary_edges %>%
  distinct(from, to, .keep_all = TRUE) %>%
  group_by(from, to) %>%
  summarise(strength = max(strength, na.rm = TRUE),
            interaction = first(interaction),
            label = first(label),
            edge_type = first(edge_type), .groups = "drop")
cat("Primary edges (measured):", nrow(primary_edges), "\n")

# 2. Cascade edges: CBL -> CIPK -> targets (literature-based)
cascade_edges <- data.frame(
  from = c("AT4G17615", "AT5G47100", "AT5G55990",  # CBLs -> CIPK23
           "AT4G17615", "AT5G47100",                 # CBLs -> CIPK9
           "AT5G55990",                               # CBL2 -> CIPK1
           "AT1G30270", "AT1G30270",                  # CIPK23 -> targets
           "AT1G01140",                               # CIPK9 -> target
           "AT3G17510",                               # CIPK1 -> target
           "AT5G37770",                               # CML24 -> SOS2 pathway
           "AT4G23650", "AT4G09570"),                 # CDPKs -> stress
  to = c("AT1G30270", "AT1G30270", "AT1G30270",
         "AT1G01140", "AT1G01140",
         "AT3G17510",
         "AT1G08090", "AT3G51860",
         "AT3G51860",
         "AT3G51860",
         "AT2G01980",
         "AT4G23650", "AT1G35720"),
  label = c(rep("phosphorylates", 3),
            rep("activates", 2),
            "activates",
            "regulates", "regulates",
            "regulates",
            "regulates",
            "activates",
            "interacts", "regulates"),
  strength = rep(NA, 13),
  interaction = rep("cascade", 13),
  edge_type = rep("cascade", 13),
  stringsAsFactors = FALSE
)

# 3. CaM/CDPK branch edges (literature-based signaling)
branch_edges <- data.frame(
  from = c("AT3G43810", "AT3G43810",   # CaM7 -> targets
           "AT3G56800",                 # CaM3 -> heat response
           "AT4G23650",                 # CDPK3 -> stress
           "AT1G35720"),                # ANNAT1 -> membrane
  to = c("AT3G51860", "AT1G08090",
         "AT4G10250",
         "AT4G10250",
         "AT4G01390"),
  label = c("regulates", "regulates",
            "induces",
            "phosphorylates",
            "interacts"),
  strength = rep(NA, 5),
  interaction = rep("branch", 5),
  edge_type = rep("branch", 5),
  stringsAsFactors = FALSE
)

# 4. LASSO gene connections
lasso_edges <- data.frame(
  from = c("Ca2+", "AT4G01390", "AT3G56800"),
  to = c("AT4G10250", "AT4G10250", "AT4G10250"),
  label = c("induces", "co-regulates", "induces"),
  strength = c(NA, NA, NA),
  interaction = c("LASSO_link", "LASSO_link", "LASSO_link"),
  edge_type = c("lasso", "lasso", "lasso"),
  stringsAsFactors = FALSE
)

# Combine all edges
edge_table <- bind_rows(primary_edges, cascade_edges, branch_edges, lasso_edges) %>%
  filter(from %in% node_defs$name & to %in% node_defs$name)
cat("Total edges:", nrow(edge_table), "\n")
cat("  Primary (measured):", sum(edge_table$edge_type == "primary"), "\n")
cat("  Cascade (literature):", sum(edge_table$edge_type == "cascade"), "\n")
cat("  Branch (literature):", sum(edge_table$edge_type == "branch"), "\n")
cat("  LASSO links:", sum(edge_table$edge_type == "lasso"), "\n")

# ============================================================
# Node table with expression
# ============================================================
expr_all <- ca2_expr %>% filter(celltype == "ALL")
expr_map <- setNames(expr_all$mean_expr, expr_all$gene)

node_table <- node_defs %>%
  mutate(expression = sapply(name, function(g) {
    if (g %in% names(expr_map)) return(expr_map[[g]])
    return(0)
  }))

cat("\nNode expression values:\n")
for (i in 1:nrow(node_table)) {
  cat("  ", node_table$alias[i], "(", node_table$name[i], "): ",
      round(node_table$expression[i], 4), "\n", sep = "")
}

# ============================================================
# Build network and apply stress-majorization layout
# ============================================================
ca2_network <- graph_from_data_frame(d = edge_table, vertices = node_table, directed = TRUE)

# Stress-majorization layout from graphlayouts package
set.seed(42)
layout_stress <- layout_with_stress(ca2_network)

# Post-layout adjustment: shift Ca2+ to top center
ca2_idx <- which(node_table$name == "Ca2+")
# Normalize layout to reasonable range
layout_stress <- layout_stress * 2
# Move Ca2+ to top
layout_stress[ca2_idx, 1] <- 0
layout_stress[ca2_idx, 2] <- max(layout_stress[, 2]) + 1.5

# Add layout to node table
node_table$x <- layout_stress[, 1]
node_table$y <- layout_stress[, 2]

# Rebuild graph with layout coordinates
ca2_network <- graph_from_data_frame(d = edge_table, vertices = node_table, directed = TRUE)

# ============================================================
# Plot with ggraph
# ============================================================
# Type colors (Phylo palette + extensions)
type_colors <- c(
  "ligand" = "#E9ED4C",
  "sensor_CBL" = "#0279EE",
  "kinase_CIPK" = "#75A025",
  "sensor_CaM" = "#FD9BED",
  "kinase_CDPK" = "#FF9400",
  "target" = "#56B4E9",
  "LASSO_gene" = "#D55E00"
)

type_labels <- c(
  "ligand" = "Ca2+ (ligand)",
  "sensor_CBL" = "CBL sensors",
  "kinase_CIPK" = "CIPK kinases",
  "sensor_CaM" = "CaM/CML sensors",
  "kinase_CDPK" = "CDPK kinases",
  "target" = "Downstream targets",
  "LASSO_gene" = "LASSO biomarkers"
)

p <- ggraph(ca2_network, layout = "manual", x = x, y = y) +
  # Primary edges (Ca2+ -> receptors, measured strength)
  geom_edge_link(aes(filter = edge_type == "primary", width = strength, alpha = strength),
                 arrow = arrow(length = unit(0.2, 'lines'), type = "closed"),
                 start_cap = circle(0.35, 'lines'), end_cap = circle(0.45, 'lines'),
                 color = "#0279EE") +
  # Cascade edges (dashed, literature)
  geom_edge_link(aes(filter = edge_type == "cascade"),
                 arrow = arrow(length = unit(0.2, 'lines'), type = "closed"),
                 start_cap = circle(0.3, 'lines'), end_cap = circle(0.4, 'lines'),
                 color = "grey50", width = 0.7, linetype = "dashed") +
  # Branch edges (dotted, literature)
  geom_edge_link(aes(filter = edge_type == "branch"),
                 arrow = arrow(length = unit(0.2, 'lines'), type = "closed"),
                 start_cap = circle(0.3, 'lines'), end_cap = circle(0.4, 'lines'),
                 color = "#FD9BED", width = 0.7, linetype = "dotted") +
  # LASSO link edges (dotted, orange)
  geom_edge_link(aes(filter = edge_type == "lasso"),
                 arrow = arrow(length = unit(0.2, 'lines'), type = "closed"),
                 start_cap = circle(0.3, 'lines'), end_cap = circle(0.4, 'lines'),
                 color = "#FF9400", width = 1, linetype = "dotted") +
  scale_edge_width(range = c(0.3, 2.5), name = "CCC strength", guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.7), name = "CCC strength") +
  # Nodes
  geom_node_point(aes(fill = type, size = expression + 0.3),
                  shape = 21, color = "grey20", stroke = 0.8) +
  geom_node_text(aes(label = alias), size = 2.8, hjust = 0.5, vjust = 0.5,
                 repel = TRUE, family = "sans", fontface = "bold",
                 bg.colour = "white", bg.r = 0.15, max.overlaps = 25) +
  scale_fill_manual(values = type_colors, name = "Node type",
                    limits = names(type_colors), labels = type_labels) +
  scale_size_continuous(range = c(3, 9), name = "Expression\n(mean)") +
  annotate(geom = "text", label = "Ca2+ Signaling Pathway (Refined)",
           x = 0, y = max(node_table$y) + 0.8, size = 5.5, fontface = "bold", family = "sans") +
  annotate(geom = "text",
           label = "Dominant CCC pathway (strength = 249.5) | 20 nodes, stress-majorization layout",
           x = 0, y = max(node_table$y) + 0.3, size = 3, family = "sans", color = "grey30") +
  coord_fixed(clip = "off") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(family = "sans", size = 8),
    legend.title = element_text(family = "sans", size = 9, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(paste0(RESULTS, "/figures/ccc_ca2_pathway_ggpathway.png"),
       plot = p, height = 9, width = 12, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ccc_ca2_pathway_ggpathway.svg"),
       plot = p, height = 9, width = 12, bg = "white")

cat("\nSaved: ccc_ca2_pathway_ggpathway.png/.svg\n")
cat("Nodes:", nrow(node_table), "| Edges:", nrow(edge_table), "\n")
cat("Layout: stress-majorization (graphlayouts::layout_with_stress)\n")
