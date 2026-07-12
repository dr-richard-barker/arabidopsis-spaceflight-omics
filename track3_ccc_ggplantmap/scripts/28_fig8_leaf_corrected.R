#!/usr/bin/env Rscript
# Redesigned Figure 8 panel a: leaf cross-sections of (i) CIPK23 and (ii) AKT1
# expression, using CORRECT gene identities (CIPK23 = AT1G30270). Replaces the
# invalid ggKEGG ath04075 panel (which highlighted CPK5 mislabeled as CIPK23).
WORK <- Sys.getenv("ASO_WORK"); if (WORK=="") WORK <- "."
suppressMessages({library(ggplot2); library(ggPlantmap); library(patchwork); library(svglite)})

roi_to_celltype <- c(
  "epidermis.stomata"="Guard","Parenchima.palisade"="Mesophyll","Parenchima.sponge"="Mesophyll",
  "epidermis.adaxial"="Epidermal","epidermis.abaxial"="Epidermal",
  "vascularbundle.xylem"="Stele","vascularbundle.phloem"="Stele","vascularbundle.bundlesheet"="Stele")

# corrected per-cell-type mean expression (from atlas re-run)
cipk23_expr <- c(Guard=0.333, Mesophyll=0.232, Epidermal=0.304, Stele=0.277)
akt1_expr   <- c(Guard=0.192, Mesophyll=0.197, Epidermal=0.367, Stele=0.288)

leaf_map <- ggPm.At.leaf.crosssection
roi_names <- unique(leaf_map$ROI.name)
mk <- function(v) data.frame(ROI.name=roi_names, val=v[roi_to_celltype[roi_names]], stringsAsFactors=FALSE)
leaf_cipk23 <- ggPlantmap.merge(leaf_map, mk(cipk23_expr), id.x="ROI.name", id.y="ROI.name")
leaf_akt1   <- ggPlantmap.merge(leaf_map, mk(akt1_expr),   id.x="ROI.name", id.y="ROI.name")

p_cipk <- ggPlantmap.heatmap(map.quant=leaf_cipk23, value.quant=val, show.legend=TRUE) +
  scale_fill_gradient(low="#ECE9E2", high="#FF9400", name="CIPK23\nexpression\n(mean)",
                      limits=c(0.20,0.36)) +
  labs(title="(i) CIPK23 expression", subtitle="Guard cells highest (0.333)") +
  theme_void(base_family="sans") +
  theme(plot.title=element_text(size=11,face="bold",hjust=0.5),
        plot.subtitle=element_text(size=9,hjust=0.5,color="#FF9400"),
        legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7))

p_akt1 <- ggPlantmap.heatmap(map.quant=leaf_akt1, value.quant=val, show.legend=TRUE) +
  scale_fill_gradient(low="#ECE9E2", high="#75A025", name="AKT1\nexpression\n(mean)",
                      limits=c(0.15,0.40)) +
  labs(title="(ii) AKT1 expression", subtitle="Epidermal highest (0.367)") +
  theme_void(base_family="sans") +
  theme(plot.title=element_text(size=11,face="bold",hjust=0.5),
        plot.subtitle=element_text(size=9,hjust=0.5,color="#75A025"),
        legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7))

p <- p_cipk + p_akt1 + plot_layout(ncol=2)
ggsave(file.path(WORK,"fig8a_leaf_corrected.png"), p, width=10, height=6, dpi=300, bg="white")
cat("Saved fig8a_leaf_corrected.png\n")
