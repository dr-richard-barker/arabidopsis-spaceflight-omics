# R dependencies for arabidopsis-spaceflight-omics
# Install with: Rscript environment.R

# CRAN packages
install.packages(c(
    "igraph",        # network analysis
    "ggraph",        # network visualization
    "ggplot2",       # plotting
    "dplyr",         # data manipulation
    "tidyr",         # data tidying
    "viridis",       # color palettes
    "cowplot",       # plot composition
    "patchwork",     # plot composition
    "gprofiler2",    # GO enrichment
    "remotes",       # GitHub installs
    "svglite",       # SVG export
    "Matrix",        # sparse matrices
    "data.table"     # fast I/O
), repos = "https://cloud.r-project.org")

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
    "AnnotationDbi",  # annotation infrastructure
    "GO.db",          # Gene Ontology database
    "clusterProfiler",# enrichment analysis (optional)
    "KEGGREST"        # KEGG pathway retrieval (Figure 8 panel a)
), update = FALSE, ask = FALSE)

# GitHub packages
remotes::install_github("leonardojo/ggPlantmap")   # spatial visualization
remotes::install_github("mrliuw/PlantCellChat")     # cell-cell communication
remotes::install_github("noriakis/ggkegg")          # KEGG pathway overlay (Figure 8 panel a)

# Package versions used in this study:
# R 4.4.x, ggplot2 4.0.3, igraph (latest), ggraph (latest),
# graphlayouts (latest), gprofiler2 (latest), ggPlantmap (GitHub HEAD),
# ggkegg (GitHub HEAD), KEGGREST (Bioconductor), png (latest)
