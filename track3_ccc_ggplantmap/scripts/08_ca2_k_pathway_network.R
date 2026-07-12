#!/usr/bin/env Rscript
# Ca2+/K+ signaling pathway network — expanded with all 6 K+ channels
# Includes AKT1 (AT2G26650, the correct gene ID), AKT2, GORK, KAT1, KAT2, KC1
# Stress-majorization layout (graphlayouts package)
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

cat("=== Ca2+/K+ Pathway Network (26 nodes, stress-majorization) ===\n\n")

# Load data
ca2_lr <- read.csv(paste0(RESULTS, "/tables/ccc_lr_pairs.csv"), stringsAsFactors = FALSE)
ca2_full <- read.csv(paste0(RESULTS, "/tables/ccc_full_communication_table.csv"), stringsAsFactors = FALSE)
ca2_strength <- ca2_full %>% filter(Signal == "Ca2+") %>%
  group_by(Interaction_name) %>% summarise(total_prob = sum(Prob, na.rm = TRUE)) %>% arrange(desc(total_prob))
k_full <- ca2_full %>% filter(Signal == "K+")
k_strength <- k_full %>% group_by(Interaction_name) %>% summarise(total_prob = sum(Prob, na.rm = TRUE)) %>% arrange(desc(total_prob))
ca2_expr <- read.csv(paste0(ATLAS, "/ca2_k_expr_per_celltype.csv"), stringsAsFactors = FALSE, check.names = FALSE)

cat("Ca2+ LR pairs:", nrow(ca2_lr), "| K+ LR pairs:", sum(ca2_lr$Signal == "K+"), "\n")

# ============================================================
# Node definitions (26 nodes)
# ============================================================
node_defs <- data.frame(
  name = c("Ca2+", "K+",
           # Ca2+ sensors (CBLs)
           "AT4G17615", "AT5G47100", "AT5G55990",
           # Ca2+ sensor kinases (CIPKs)
           "AT1G30270", "AT3G17510", "AT1G01140",
           # CaM/CML sensors
           "AT1G76650", "AT4G20780", "AT5G37770", "AT3G43810", "AT3G56800",
           # CDPKs
           "AT4G23650", "AT4G09570",
           # Ca2+ downstream targets
           "AT3G51860", "AT1G08090", "AT2G01980", "AT1G35720",
           # K+ channels (ALL 6)
           "AT2G26650", "AT4G22200", "AT5G37500", "AT5G46240", "AT4G18290", "AT4G32650",
           # LASSO genes
           "AT4G10250", "AT4G01390"),
  alias = c("Ca2+", "K+",
            "CBL1", "CBL9", "CBL2",
            "CIPK23", "CIPK1", "CIPK9",
            "CML38", "CML42", "CML24", "CaM7", "CaM3",
            "CDPK3", "CDPK4",
            "CAX3", "NRT2.1", "SOS1", "ANNAT1",
            "AKT1", "AKT2", "GORK", "KAT1", "KAT2", "KC1",
            "HSP22.0\n(LASSO+)", "TRAF-like\n(LASSO-)"),
  type = c("ligand_ca2", "ligand_k",
           rep("sensor_CBL", 3),
           rep("kinase_CIPK", 3),
           rep("sensor_CaM", 5),
           rep("kinase_CDPK", 2),
           rep("target_ca2", 4),
           rep("channel_K", 6),
           rep("LASSO_gene", 2)),
  stringsAsFactors = FALSE
)
cat("Nodes:", nrow(node_defs), "\n")

# ============================================================
# Edge table
# ============================================================
node_genes <- setdiff(node_defs$name, c("Ca2+", "K+"))

# 1. Primary Ca2+ edges (measured CCC)
primary_ca2 <- data.frame()
ca2_lr_db <- ca2_lr %>% filter(Signal == "Ca2+")
for (i in 1:nrow(ca2_lr_db)) {
  rec <- ca2_lr_db$Receptor[i]; lig <- ca2_lr_db$Ligand[i]
  interaction_name <- ca2_lr_db$Interaction_name[i]
  strength_val <- ca2_strength$total_prob[ca2_strength$Interaction_name == interaction_name]
  if (length(strength_val) == 0) next
  if (rec %in% node_genes) {
    primary_ca2 <- rbind(primary_ca2, data.frame(
      from="Ca2+", to=rec, label="binds", strength=strength_val,
      interaction=interaction_name, edge_type="primary_ca2", pathway="Ca2+", stringsAsFactors=FALSE))
  }
  if (lig %in% node_genes && rec %in% node_genes && lig != rec) {
    primary_ca2 <- rbind(primary_ca2, data.frame(
      from=lig, to=rec, label="interacts", strength=strength_val,
      interaction=interaction_name, edge_type="primary_ca2", pathway="Ca2+", stringsAsFactors=FALSE))
  }
}
primary_ca2 <- primary_ca2 %>% distinct(from, to, .keep_all=TRUE) %>%
  group_by(from, to) %>% summarise(strength=max(strength, na.rm=TRUE),
    interaction=first(interaction), label=first(label), edge_type=first(edge_type), pathway=first(pathway), .groups="drop")
cat("Primary Ca2+ edges:", nrow(primary_ca2), "\n")

# 2. Primary K+ edges (measured CCC)
primary_k <- data.frame()
k_lr_db <- ca2_lr %>% filter(Signal == "K+")
for (i in 1:nrow(k_lr_db)) {
  rec <- k_lr_db$Receptor[i]; lig <- k_lr_db$Ligand[i]
  interaction_name <- k_lr_db$Interaction_name[i]
  strength_val <- k_strength$total_prob[k_strength$Interaction_name == interaction_name]
  if (length(strength_val) == 0) next
  if (rec %in% node_genes) {
    primary_k <- rbind(primary_k, data.frame(
      from="K+", to=rec, label="binds", strength=strength_val,
      interaction=interaction_name, edge_type="primary_k", pathway="K+", stringsAsFactors=FALSE))
  }
  if (lig %in% node_genes && rec %in% node_genes && lig != rec) {
    primary_k <- rbind(primary_k, data.frame(
      from=lig, to=rec, label="interacts", strength=strength_val,
      interaction=interaction_name, edge_type="primary_k", pathway="K+", stringsAsFactors=FALSE))
  }
}
primary_k <- primary_k %>% distinct(from, to, .keep_all=TRUE) %>%
  group_by(from, to) %>% summarise(strength=max(strength, na.rm=TRUE),
    interaction=first(interaction), label=first(label), edge_type=first(edge_type), pathway=first(pathway), .groups="drop")
cat("Primary K+ edges:", nrow(primary_k), "\n")

# 3. Cascade edges (literature: CBL->CIPK->targets)
cascade_edges <- data.frame(
  from = c("AT4G17615","AT5G47100","AT5G55990",  # CBLs -> CIPK23
           "AT4G17615","AT5G47100",               # CBLs -> CIPK9
           "AT5G55990",                            # CBL2 -> CIPK1
           "AT1G30270",                            # CIPK23 -> AKT1 (KEY crosstalk)
           "AT1G30270",                            # CIPK23 -> NRT2.1
           "AT1G30270",                            # CIPK23 -> CAX3
           "AT1G01140",                            # CIPK9 -> CAX3
           "AT3G17510",                            # CIPK1 -> CAX3
           "AT5G37770",                            # CML24 -> SOS1
           "AT4G23650","AT4G09570"),               # CDPKs -> targets
  to = c("AT1G30270","AT1G30270","AT1G30270",
         "AT1G01140","AT1G01140",
         "AT3G17510",
         "AT2G26650",
         "AT1G08090",
         "AT3G51860",
         "AT3G51860",
         "AT3G51860",
         "AT2G01980",
         "AT4G23650","AT1G35720"),
  label = c(rep("phosphorylates",3), rep("activates",2), "activates",
            "phosphorylates", "regulates", "regulates", "regulates", "regulates",
            "activates", "interacts", "regulates"),
  strength = rep(NA, 14), interaction = rep("cascade", 14),
  edge_type = rep("cascade", 14), pathway = rep("Ca2+", 14),
  stringsAsFactors = FALSE
)

# 4. CaM/CDPK branch edges
branch_edges <- data.frame(
  from = c("AT3G43810","AT3G43810",   # CaM7 -> targets
           "AT3G56800",                 # CaM3 -> heat response
           "AT4G23650",                 # CDPK3 -> stress
           "AT1G35720"),                # ANNAT1 -> membrane
  to = c("AT3G51860","AT1G08090",
         "AT4G10250","AT4G10250","AT4G01390"),
  label = c("regulates","regulates","induces","phosphorylates","interacts"),
  strength = rep(NA, 5), interaction = rep("branch", 5),
  edge_type = rep("branch", 5), pathway = rep("Ca2+", 5),
  stringsAsFactors = FALSE
)

# 5. K+ channel interaction edges (literature + measured)
k_channel_edges <- data.frame(
  from = c("AT4G32650",   # KC1 -> AKT1 (measured K+ CCC, strongest)
           "AT4G22200",   # AKT2 -> AKT1 (measured K+ CCC)
           "AT5G47100",   # CBL9 -> AKT1 (CROSSTALK: Ca2+ sensor regulates K+ channel)
           "AT4G18290",   # KAT2 -> AKT1 (measured K+ CCC)
           "AT5G37500"),  # GORK -> AKT2 (measured K+ CCC)
  to = c("AT2G26650","AT2G26650","AT2G26650","AT2G26650","AT4G22200"),
  label = c("interacts","interacts","regulates","interacts","interacts"),
  strength = c(5.01, 1.41, NA, 2.16, 1.01),
  interaction = c("K_CCC","K_CCC","crosstalk","K_CCC","K_CCC"),
  edge_type = c("k_channel","k_channel","crosstalk","k_channel","k_channel"),
  pathway = c("K+","K+","Ca2+/K+","K+","K+"),
  stringsAsFactors = FALSE
)

# 6. LASSO gene connections
lasso_edges <- data.frame(
  from = c("Ca2+","AT4G01390","AT3G56800"),
  to = c("AT4G10250","AT4G10250","AT4G10250"),
  label = c("induces","co-regulates","induces"),
  strength = c(NA,NA,NA), interaction = rep("LASSO_link",3),
  edge_type = rep("lasso",3), pathway = rep("Ca2+",3),
  stringsAsFactors = FALSE
)

# Combine all edges
edge_table <- bind_rows(primary_ca2, primary_k, cascade_edges, branch_edges, k_channel_edges, lasso_edges) %>%
  filter(from %in% node_defs$name & to %in% node_defs$name)
cat("Total edges:", nrow(edge_table), "\n")
cat("  Ca2+ primary:", sum(edge_table$edge_type=="primary_ca2"), "\n")
cat("  K+ primary:", sum(edge_table$edge_type=="primary_k"), "\n")
cat("  Cascade:", sum(edge_table$edge_type=="cascade"), "\n")
cat("  Branch:", sum(edge_table$edge_type=="branch"), "\n")
cat("  K+ channel:", sum(edge_table$edge_type=="k_channel"), "\n")
cat("  Crosstalk:", sum(edge_table$edge_type=="crosstalk"), "\n")
cat("  LASSO:", sum(edge_table$edge_type=="lasso"), "\n")

# ============================================================
# Node table with expression
# ============================================================
expr_all <- ca2_expr %>% filter(celltype == "ALL")
expr_map <- setNames(expr_all$mean_expr, expr_all$gene)
node_table <- node_defs %>% mutate(expression = sapply(name, function(g) {
  if (g %in% names(expr_map)) return(expr_map[[g]]); return(0) }))

# ============================================================
# Build network with stress-majorization layout
# ============================================================
ca2_k_network <- graph_from_data_frame(d = edge_table, vertices = node_table, directed = TRUE)
set.seed(42)
layout_stress <- layout_with_stress(ca2_k_network)
layout_stress <- layout_stress * 2

# Shift Ca2+ and K+ to top
ca2_idx <- which(node_table$name == "Ca2+")
k_idx <- which(node_table$name == "K+")
layout_stress[ca2_idx, 1] <- -1.5
layout_stress[ca2_idx, 2] <- max(layout_stress[, 2]) + 1.5
layout_stress[k_idx, 1] <- 1.5
layout_stress[k_idx, 2] <- max(layout_stress[, 2]) + 1.5

node_table$x <- layout_stress[, 1]
node_table$y <- layout_stress[, 2]
ca2_k_network <- graph_from_data_frame(d = edge_table, vertices = node_table, directed = TRUE)

# ============================================================
# Plot
# ============================================================
type_colors <- c(
  "ligand_ca2" = "#E9ED4C", "ligand_k" = "#75A025",
  "sensor_CBL" = "#0279EE", "kinase_CIPK" = "#0279EE",
  "sensor_CaM" = "#FD9BED", "kinase_CDPK" = "#FF9400",
  "target_ca2" = "#56B4E9", "channel_K" = "#009E73",
  "LASSO_gene" = "#D55E00"
)
type_labels <- c(
  "ligand_ca2" = "Ca2+ (ligand)", "ligand_k" = "K+ (ligand)",
  "sensor_CBL" = "CBL sensors", "kinase_CIPK" = "CIPK kinases",
  "sensor_CaM" = "CaM/CML sensors", "kinase_CDPK" = "CDPK kinases",
  "target_ca2" = "Ca2+ targets", "channel_K" = "K+ channels",
  "LASSO_gene" = "LASSO biomarkers"
)

p <- ggraph(ca2_k_network, layout = "manual", x = x, y = y) +
  # Ca2+ primary edges (blue)
  geom_edge_link(aes(filter = edge_type == "primary_ca2", width = strength, alpha = strength),
    arrow = arrow(length = unit(0.18, 'lines'), type = "closed"),
    start_cap = circle(0.3, 'lines'), end_cap = circle(0.4, 'lines'), color = "#0279EE") +
  # K+ primary edges (green)
  geom_edge_link(aes(filter = edge_type == "primary_k", width = strength, alpha = strength),
    arrow = arrow(length = unit(0.18, 'lines'), type = "closed"),
    start_cap = circle(0.3, 'lines'), end_cap = circle(0.4, 'lines'), color = "#009E73") +
  # Cascade edges (dashed grey)
  geom_edge_link(aes(filter = edge_type == "cascade"),
    arrow = arrow(length = unit(0.18, 'lines'), type = "closed"),
    start_cap = circle(0.25, 'lines'), end_cap = circle(0.35, 'lines'),
    color = "grey50", width = 0.6, linetype = "dashed") +
  # Branch edges (dotted pink)
  geom_edge_link(aes(filter = edge_type == "branch"),
    arrow = arrow(length = unit(0.18, 'lines'), type = "closed"),
    start_cap = circle(0.25, 'lines'), end_cap = circle(0.35, 'lines'),
    color = "#FD9BED", width = 0.6, linetype = "dotted") +
  # K+ channel edges (dashed green)
  geom_edge_link(aes(filter = edge_type == "k_channel", width = strength, alpha = strength),
    arrow = arrow(length = unit(0.18, 'lines'), type = "closed"),
    start_cap = circle(0.25, 'lines'), end_cap = circle(0.35, 'lines'),
    color = "#009E73", width = 0.8, linetype = "dashed") +
  # Crosstalk edges (orange, bold)
  geom_edge_link(aes(filter = edge_type == "crosstalk"),
    arrow = arrow(length = unit(0.2, 'lines'), type = "closed"),
    start_cap = circle(0.25, 'lines'), end_cap = circle(0.35, 'lines'),
    color = "#FF9400", width = 1.5, linetype = "longdash") +
  # LASSO edges (dotted orange)
  geom_edge_link(aes(filter = edge_type == "lasso"),
    arrow = arrow(length = unit(0.18, 'lines'), type = "closed"),
    start_cap = circle(0.25, 'lines'), end_cap = circle(0.35, 'lines'),
    color = "#D55E00", width = 0.8, linetype = "dotted") +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.7), guide = "none") +
  geom_node_point(aes(fill = type, size = expression + 0.3), shape = 21, color = "grey20", stroke = 0.8) +
  geom_node_text(aes(label = alias), size = 2.5, hjust = 0.5, vjust = 0.5,
    repel = TRUE, family = "sans", fontface = "bold", bg.colour = "white", bg.r = 0.15, max.overlaps = 30) +
  scale_fill_manual(values = type_colors, name = "Node type", limits = names(type_colors), labels = type_labels) +
  scale_size_continuous(range = c(3, 9), name = "Expression\n(mean)") +
  annotate(geom = "text", label = "Ca2+/K+ Signaling Network",
    x = 0, y = max(node_table$y) + 0.8, size = 5.5, fontface = "bold", family = "sans") +
  annotate(geom = "text",
    label = "26 nodes | Ca2+ (blue) + K+ (green) + crosstalk (orange) | stress-majorization layout",
    x = 0, y = max(node_table$y) + 0.3, size = 3, family = "sans", color = "grey30") +
  coord_fixed(clip = "off") + theme_void() +
  theme(legend.position = "right", legend.text = element_text(family = "sans", size = 8),
    legend.title = element_text(family = "sans", size = 9, face = "bold"), plot.margin = margin(20, 20, 20, 20))

ggsave(paste0(RESULTS, "/figures/ccc_ca2_k_pathway_ggpathway.png"), p, height = 10, width = 13, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/ccc_ca2_k_pathway_ggpathway.svg"), p, height = 10, width = 13, bg = "white")
cat("\nSaved: ccc_ca2_k_pathway_ggpathway.png/.svg\n")
cat("Nodes:", nrow(node_table), "| Edges:", nrow(edge_table), "\n")
