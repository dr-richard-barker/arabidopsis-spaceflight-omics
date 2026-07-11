#!/usr/bin/env Rscript
# GO enrichment analysis for 85 LASSO panel genes using g:Profiler
# Custom 2000-gene background, BH-FDR correction
# Three queries: all 85, 41 flight-up, 44 ground-up
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
library(gprofiler2)
library(ggplot2)
library(dplyr)
library(tidyr)

cat("=== GO Enrichment Analysis (85 LASSO Genes) ===\n\n")

# Load data
biomarker <- read.csv(paste0(RESULTS, "/tables/biomarker_panel.csv"), stringsAsFactors = FALSE)
background <- read.csv(paste0(RESULTS, "/tables/all_feature_stability.csv"), stringsAsFactors = FALSE)$feature

genes_all <- biomarker$feature
genes_pos <- biomarker$feature[biomarker$mean_coefficient > 0]  # flight-up
genes_neg <- biomarker$feature[biomarker$mean_coefficient < 0]  # ground-up

cat("Query sizes: all=", length(genes_all), ", flight-up=", length(genes_pos),
    ", ground-up=", length(genes_neg), "\n")
cat("Background size:", length(background), "\n\n")

# Run g:Profiler for all three queries
sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG")

run_enrichment <- function(genes, query_name, bg) {
  cat("Running enrichment for:", query_name, "(", length(genes), "genes )\n")
  res <- gost(query = genes, organism = "athaliana", sources = sources,
              significant = FALSE, correction_method = "fdr",
              custom_bg = bg, evcodes = TRUE)
  if (is.null(res$result) || nrow(res$result) == 0) {
    cat("  No results returned\n")
    return(NULL)
  }
  res$result$query_name <- query_name
  res$result
}

res_all <- run_enrichment(genes_all, "All_85", background)
res_pos <- run_enrichment(genes_pos, "Flight_up_41", background)
res_neg <- run_enrichment(genes_neg, "Ground_up_44", background)

# Combine all results
all_results <- bind_rows(res_all, res_pos, res_neg)
cat("\nTotal terms tested across all queries:", nrow(all_results), "\n")

# Save full results
full_out <- all_results %>%
  select(query_name, source, term_id, term_name, p_value, term_size,
         query_size, intersection_size, precision, recall, intersection) %>%
  arrange(query_name, p_value)
write.csv(full_out, paste0(RESULTS, "/tables/go_enrichment_85_genes.csv"), row.names = FALSE)
cat("Saved full results: go_enrichment_85_genes.csv\n")

# Significant terms (padj < 0.05)
sig_results <- all_results %>% filter(p_value < 0.05) %>% arrange(query_name, p_value)
cat("\nSignificant terms (BH-FDR < 0.05):\n")
print(sig_results[, c("query_name", "source", "term_name", "p_value", "intersection_size")])
write.csv(sig_results[, c("query_name", "source", "term_id", "term_name", "p_value",
                           "term_size", "query_size", "intersection_size", "intersection")],
          paste0(RESULTS, "/tables/go_enrichment_significant.csv"), row.names = FALSE)
cat("Saved significant results: go_enrichment_significant.csv\n")

# Create summary table
summary_table <- sig_results %>%
  mutate(neg_log10_p = -log10(p_value)) %>%
  select(query_name, source, term_name, p_value, neg_log10_p,
         intersection_size, term_size, intersection) %>%
  arrange(desc(neg_log10_p))
write.csv(summary_table, paste0(RESULTS, "/tables/go_enrichment_summary.csv"), row.names = FALSE)

# ============================================================
# Visualization 1: Dot plot of top significant terms
# ============================================================
cat("\n=== Creating dot plot ===\n")

# Prepare data for dot plot - top 15 terms by significance across all queries
dot_data <- sig_results %>%
  mutate(neg_log10_p = -log10(p_value)) %>%
  mutate(query_label = case_when(
    query_name == "All_85" ~ "All 85 genes",
    query_name == "Flight_up_41" ~ "Flight-up (41)",
    query_name == "Ground_up_44" ~ "Ground-up (44)"
  )) %>%
  mutate(term_short = ifelse(nchar(term_name) > 50, paste0(substr(term_name, 1, 47), "..."), term_name))

# Get top 15 unique terms by min p-value
top_terms <- dot_data %>%
  group_by(source, term_name) %>%
  summarise(min_p = min(p_value), .groups = "drop") %>%
  arrange(min_p) %>%
  head(15)

dot_data_top <- dot_data %>%
  inner_join(top_terms %>% select(source, term_name), by = c("source", "term_name")) %>%
  mutate(term_short = factor(term_short, levels = rev(unique(term_short[order(term_short)]))))

# Also include non-significant terms for the queries where they were tested (show as small dots)
all_dot <- all_results %>%
  filter(source %in% c("GO:BP", "GO:MF", "GO:CC", "KEGG")) %>%
  mutate(neg_log10_p = -log10(p_value)) %>%
  mutate(query_label = case_when(
    query_name == "All_85" ~ "All 85 genes",
    query_name == "Flight_up_41" ~ "Flight-up (41)",
    query_name == "Ground_up_44" ~ "Ground-up (44)"
  )) %>%
  inner_join(top_terms %>% select(source, term_name), by = c("source", "term_name")) %>%
  mutate(term_short = ifelse(nchar(term_name) > 50, paste0(substr(term_name, 1, 47), "..."), term_name)) %>%
  mutate(term_short = factor(term_short, levels = rev(unique(top_terms$term_name[order(top_terms$min_p)])))) %>%
  mutate(significant = p_value < 0.05)

p_dot <- ggplot(all_dot, aes(x = query_label, y = term_short)) +
  geom_point(aes(size = intersection_size, fill = neg_log10_p), shape = 21, color = "grey30", stroke = 0.5) +
  scale_fill_gradient(low = "#FFF5EE", high = "#D55E00", name = "-log10(padj)") +
  scale_size_continuous(range = c(2, 8), name = "Gene count") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title = "GO/KEGG Enrichment of LASSO Biomarker Panel",
       subtitle = "Custom 2000-gene background, BH-FDR correction") +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.grid.major = element_line(color = "grey92"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(paste0(RESULTS, "/figures/go_enrichment_dotplot.png"), p_dot,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/go_enrichment_dotplot.svg"), p_dot,
       width = 10, height = 8, bg = "white")
cat("Saved: go_enrichment_dotplot.png/.svg\n")

# ============================================================
# Visualization 2: Directed comparison (flight-up vs ground-up)
# ============================================================
cat("\n=== Creating directed comparison bar plot ===\n")

# Get significant terms for each direction
sig_pos <- sig_results %>% filter(query_name == "Flight_up_41") %>%
  mutate(direction = "Flight-up") %>% head(10)
sig_neg <- sig_results %>% filter(query_name == "Ground_up_44") %>%
  mutate(direction = "Ground-up") %>% head(10)

# Also include near-significant ground-up terms (p < 0.1) to show the contrast
near_neg <- all_results %>% filter(query_name == "Ground_up_44", p_value < 0.1) %>%
  mutate(direction = "Ground-up") %>% head(5)

directed_data <- bind_rows(
  sig_pos %>% mutate(direction = "Flight-up"),
  sig_neg %>% mutate(direction = "Ground-up")
) %>%
  mutate(neg_log10_p = -log10(p_value)) %>%
  mutate(term_short = ifelse(nchar(term_name) > 45, paste0(substr(term_name, 1, 42), "..."), term_name)) %>%
  arrange(direction, desc(neg_log10_p))

# Create the bar plot
p_directed <- ggplot(directed_data, aes(x = neg_log10_p, y = reorder(term_short, neg_log10_p),
                                         fill = direction)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = paste0(intersection_size, " genes")), hjust = -0.1, size = 3, family = "sans") +
  scale_fill_manual(values = c("Flight-up" = "#D55E00", "Ground-up" = "#0072B2"),
                    name = "Direction") +
  facet_wrap(~ direction, scales = "free_y", ncol = 2) +
  labs(x = "-log10(adjusted p-value)", y = NULL,
       title = "Directed GO/KEGG Enrichment: Flight-up vs Ground-up Genes",
       subtitle = "Flight-up: stress/heat/ER-proteostasis | Ground-up: flavonoid metabolism") +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  ) +
  xlim(0, max(directed_data$neg_log10_p) * 1.3)

ggsave(paste0(RESULTS, "/figures/go_enrichment_directed.png"), p_directed,
       width = 14, height = 6, dpi = 300, bg = "white")
ggsave(paste0(RESULTS, "/figures/go_enrichment_directed.svg"), p_directed,
       width = 14, height = 6, bg = "white")
cat("Saved: go_enrichment_directed.png/.svg\n")

# Print final summary
cat("\n=== ENRICHMENT SUMMARY ===\n")
cat("All 85 genes:", sum(sig_results$query_name == "All_85"), "significant terms\n")
cat("Flight-up (41):", sum(sig_results$query_name == "Flight_up_41"), "significant terms\n")
cat("Ground-up (44):", sum(sig_results$query_name == "Ground_up_44"), "significant terms\n")
cat("\nDone.\n")
