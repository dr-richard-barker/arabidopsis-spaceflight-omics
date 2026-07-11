#!/usr/bin/env Rscript
# Export metadata and gene names as CSV (complement to the HDF5 sparse matrix)
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

library(Seurat)

seu <- readRDS(paste0(ATLAS, "/seedling_6d.rds"))
DefaultAssay(seu) <- "RNA"
counts <- GetAssayData(seu, layer = "counts")

# Export gene names
genes <- rownames(counts)
write.csv(data.frame(gene = genes), paste0(ATLAS, "/gene_names.csv"), row.names = FALSE)

# Export cell metadata
meta <- seu@meta.data
meta$cell_id <- rownames(meta)
write.csv(meta, paste0(ATLAS, "/cell_metadata.csv"), row.names = FALSE)

# Export UMAP
if ("umap" %in% names(seu@reductions)) {
    umap <- Embeddings(seu, reduction = "umap")
    umap_df <- data.frame(cell_id = rownames(umap), UMAP_1 = umap[,1], UMAP_2 = umap[,2])
    write.csv(umap_df, paste0(ATLAS, "/umap_coords.csv"), row.names = FALSE)
}

cat("✓ Metadata exported\n")
cat("  Genes:", length(genes), "\n")
cat("  Cells:", nrow(meta), "\n")
cat("  CellTypes:", paste(unique(meta$CellType), collapse = ", "), "\n")
