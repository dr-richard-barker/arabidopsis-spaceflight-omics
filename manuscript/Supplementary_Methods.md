# Supplementary Methods

## Spaceflight Biomarker Discovery and Single-Cell Atlas Integration in *Arabidopsis thaliana*

---

## S1. OSDR Data Retrieval and Harmonization

### S1.1 Data Source
Gene expression data and metadata were retrieved from the NASA Open Science Data Repository (OSDR; https://visualization.osdr.nasa.gov/biodata/api/). Six *Arabidopsis thaliana* spaceflight studies were selected based on the following criteria: (1) whole-genome transcriptomic data (microarray or RNA-seq), (2) flight vs. ground control comparison, (3) vegetative tissue (seedling, leaf, root, or whole plant), and (4) raw or processed expression data available via the OSDR API.

### S1.2 Studies Included
| OSDR ID | Mission | Tissue | Platform | Samples |
|---------|---------|--------|----------|--------|
| OSD-37 | ISS | Seedling | Microarray | 24 |
| OSD-678 | ISS | Leaf | RNA-seq | 24 |
| OSD-38 | ISS | Root | Microarray | 24 |
| OSD-321 | ISS | Seedling | Microarray | 24 |
| OSD-120 | Shuttle | Seedling | Microarray | 24 |
| OSD-624 | ISS | Whole plant | RNA-seq | 36 |

Total: 156 samples (78 flight, 78 ground).

### S1.3 Data Processing
- **Microarray data**: Processed expression matrices were retrieved directly from OSDR. Probes were mapped to AGI gene identifiers (ATXGXXXXX) using the TAIR10 annotation. Probes mapping to multiple genes were excluded; probes mapping to the same gene were collapsed by taking the maximum expression value.
- **RNA-seq data**: Count matrices were retrieved from OSDR. Genes with <10 counts across all samples were filtered. Counts were normalized using DESeq2's median-of-ratios method and variance-stabilizing transformation (VST) was applied for downstream analysis.
- **Harmonization**: All expression matrices were converted to log2-transformed values. Gene identifiers were harmonized to AGI format. Samples were annotated with condition (flight/ground), study ID, tissue type, and mission. The final pooled matrix contained 32,548 unique genes across 156 samples.

### S1.4 Code
```python
# OSDR data retrieval (simplified)
import requests, pandas as pd

def retrieve_osdr_study(study_id):
    url = f"https://visualization.osdr.nasa.gov/biodata/api/dataset/{study_id}"
    response = requests.get(url)
    # Parse ISA-Tab metadata and expression data
    # ... (full implementation in track1_lasso_osdr/scripts/01_retrieve_osdr_data.py)
    return expression_matrix, sample_metadata

studies = ["OSD-37", "OSD-678", "OSD-38", "OSD-321", "OSD-120", "OSD-624"]
all_data = []
for sid in studies:
    expr, meta = retrieve_osdr_study(sid)
    all_data.append((sid, expr, meta))
```

---

## S2. LASSO Biomarker Panel Construction

### S2.1 Feature Pre-selection
The top 2,000 most variable genes were selected from the pooled 32,548-gene matrix using variance stabilizing transformation followed by selection of the top features by dispersion. This pre-selection step reduces the feature space to a manageable size for regularized regression while retaining genes with biological variability across conditions.

### S2.2 LASSO Regression with Study-Stratified Cross-Validation
A LASSO (Least Absolute Shrinkage and Selection Operator) logistic regression model was trained to predict flight vs. ground condition. To account for study-specific batch effects and assess generalizability, we employed a **leave-one-study-out cross-validation (LOOS-CV)** strategy:

1. For each of the 6 studies, the study was held out as the test set.
2. The remaining 5 studies were used for training, with an inner 5-fold CV to select the optimal LASSO regularization parameter (λ).
3. The model was evaluated on the held-out study, recording AUC, sensitivity, and specificity.
4. This process was repeated 30 times with different random seeds for the inner CV to assess stability.

### S2.3 Stability Selection
Features selected in ≥50% of the 30×6 = 180 CV iterations were defined as "stable" biomarkers. The selection frequency and mean coefficient (with standard deviation) were recorded for each gene. This yielded 85 stable features constituting the final biomarker panel.

### S2.4 Key Results
- Mean CV AUC: 0.734 (95% CI: 0.431–1.0)
- Mean sensitivity: 0.845, Mean specificity: 0.670
- 85 stable features at ≥50% selection frequency
- Top 4 genes (100% selection): AT4G10250 (HSP22.0, coef=+0.723), AT3G07365 (NAT, coef=+0.828), AT2G14247 (chloroplast protein, coef=−0.364), AT4G01390 (TRAF-like, coef=−0.819)

### S2.5 Code
```r
# LASSO with LOOS-CV (simplified)
library(glmnet)
studies <- unique(metadata$study)
all_coefs <- list()
for (held_out in studies) {
  train_idx <- metadata$study != held_out
  test_idx <- metadata$study == held_out
  cv_fit <- cv.glmnet(X[train_idx,], y[train_idx], family="binomial", alpha=1)
  pred <- predict(cv_fit, X[test_idx,], s="lambda.min", type="response")
  auc <- pROC::auc(y[test_idx], pred)
  coefs <- coef(cv_fit, s="lambda.min")
  all_coefs[[held_out]] <- coefs
}
# Stability: count selection frequency across all folds
```

---

## S3. Atlas Conversion (Seurat RDS to AnnData)

### S3.1 Data Source
The *Arabidopsis* seedling single-cell atlas was obtained from Lee, Nobori, Illouz-Eliaz et al. (Nature Plants, 2025; GEO accession GSE226097). The atlas contains 41,314 cells across five major cell types: Epidermal (19,233), Mesophyll (14,952), Stele (4,240), Meristematic (2,800), and Guard (89).

### S3.2 Conversion Process
The atlas was provided as a Seurat RDS object. Conversion to AnnData format (for Python/scanpy compatibility) was performed using the `sceasy` and `anndata` packages:

1. The Seurat object was loaded in R.
2. The count matrix, cell metadata (including `CellType` annotations), and gene metadata were extracted.
3. UMAP coordinates were extracted from the Seurat reduction slot.
4. Data was exported to HDF5 format and read into Python using `anndata.read_h5ad()`.
5. Quality control: cells with <200 genes and genes expressed in <3 cells were filtered.
6. Normalization: `sc.pp.normalize_total(adata, target_sum=1e4)` followed by `sc.pp.log1p(adata)`.
7. Highly variable genes: `sc.pp.highly_variable_genes(adata, n_top_genes=3000)`.
8. The final AnnData object contained 41,314 cells × 22,375 genes.

### S3.3 Subsampling
For computational efficiency in downstream scPlantLLM inference, the atlas was subsampled to 9,998 cells using stratified sampling proportional to cell type frequencies.

### S3.4 Code
```r
# Seurat to AnnData conversion (simplified)
library(Seurat)
library(sceasy)
seurat_obj <- readRDS("seedling_6d.rds")
# Export counts, metadata, and embeddings
# ... (full implementation in track2_scplantllm_atlas/scripts/01_convert_atlas.R)
```

```python
# Load and process in Python
import scanpy as sc
adata = sc.read_h5ad("seedling_6d_anndata.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
```

---

## S4. scPlantLLM Zero-Shot Inference

### S4.1 Model Description
scPlantLLM (Cao et al., Genomics Proteomics Bioinformatics, 2025) is a foundation model for plant single-cell RNA-seq data, pre-trained on ~1.7 million plant cells across multiple species and tissues. The model uses a Transformer encoder architecture with a masked language modeling objective, taking binned gene expression values as input tokens.

### S4.2 Pretrained Weight Conversion (flash_attn to PyTorch)

The pretrained scPlantLLM weights were stored in the `flash_attn` FlashMHA format, which uses a combined `Wqkv` (query-key-value) projection. The standard PyTorch `nn.TransformerEncoderLayer` expects separate `in_proj_weight` and `in_proj_bias` parameters. A conversion was required to load the weights into a standard PyTorch model.

**Key differences:**
- `flash_attn` FlashMHA: `self_attn.Wqkv.weight` [1536, 512] and `self_attn.Wqkv.bias` [1536]
- PyTorch TransformerEncoderLayer: `self_attn.in_proj_weight` [1536, 512] and `self_attn.in_proj_bias` [1536]

The layout is identical: both concatenate Q, K, V projections along the first dimension (512×3 = 1536). The conversion is a direct key rename.

**Additional modifications:**
1. `grad_reverse_discriminator.*` keys (10 keys, used during pre-training for domain adaptation) were skipped — these are not present in the inference-only model.
2. The `value_encoder.embedding.weight` has shape [103, 512], requiring `n_input_bins=103` (not 101) when constructing the model. The extra 2 bins account for the CLS token and padding token.

**Conversion code:**
```python
import torch

# Load original flash_attn weights
sd = torch.load('model_params/scPlantLLM_model.pth', map_location='cpu')
print(f"Original keys: {len(sd)}")  # 134 keys

new_sd = {}
for k, v in sd.items():
    # Convert Wqkv.weight -> in_proj_weight (same [1536, 512] layout)
    if 'Wqkv.weight' in k:
        new_key = k.replace('Wqkv.weight', 'in_proj_weight')
        new_sd[new_key] = v  # shape: [1536, 512] = [3*512, 512] (Q,K,V concatenated)
    # Convert Wqkv.bias -> in_proj_bias
    elif 'Wqkv.bias' in k:
        new_key = k.replace('Wqkv.bias', 'in_proj_bias')
        new_sd[new_key] = v  # shape: [1536]
    # Skip gradient reversal discriminator (pre-training only)
    elif 'grad_reverse_discriminator' in k:
        continue
    # Keep all other parameters unchanged
    else:
        new_sd[k] = v

torch.save(new_sd, 'model_params/scPlantLLM_model_converted.pth')
print(f"Converted keys: {len(new_sd)}")  # 124 keys
```

**Verification:**
- Original: 134 keys (12 Wqkv + 10 grad_reverse_discriminator + 112 others)
- Converted: 124 keys (12 in_proj + 112 others)
- 0 Wqkv keys remaining, 0 grad_reverse_discriminator keys remaining
- `value_encoder.embedding.weight`: [103, 512] (requires `n_input_bins=103`)
- `encoder.embedding.weight`: [185622, 512] (gene vocabulary size)
- All 124 parameters loaded successfully into the model with `strict=True`

### S4.3 Model Architecture and Inference Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| ntoken | 185,622 | Gene vocabulary size |
| d_model | 512 | Hidden dimension |
| nhead | 8 | Number of attention heads |
| d_hid | 512 | Feedforward hidden dimension |
| nlayers | 6 | Transformer encoder layers |
| nlayers_cls | 3 | Classification head layers |
| n_cls | 44 | Number of cell type classes |
| n_input_bins | 103 | Input bins (101 expression + CLS + padding) |
| CLS_ID | 185,621 | CLS token ID |
| PAD_VALUE | −2.0 | Padding value for expression |
| PAD_TOKEN_ID | 0 | Padding token ID |
| N_BINS | 101 | Number of expression bins |
| seq_length | 500 | Sequence length (reduced from 1500 for CPU) |
| batch_size | 16 | Batch size |

### S4.4 Inference Procedure

1. **Gene filtering**: Only genes present in both the atlas and the scPlantLLM vocabulary (97.4% of atlas genes) were retained.
2. **Binning**: Expression values were binned into 101 bins using the model's binning function.
3. **Sequence construction**: For each cell, the top 500 most highly expressed genes were selected and arranged as a sequence of (gene_id, bin_value) pairs. A CLS token was prepended.
4. **Category mode**: In category mode, padding values (−2.0) were replaced with `N_BINS` (101) before the forward pass.
5. **Forward pass**: The sequence was passed through the Transformer encoder. The CLS token embedding was extracted and passed through the classification head to obtain cell type probabilities.
6. **Zero-shot prediction**: The predicted cell type was the argmax of the output probabilities. No fine-tuning was performed.

### S4.5 Key Results
- Exact accuracy: 19.0%
- Loose mapping accuracy: 28.1% (2810/9998 cells)
- Per-type (loose): Mesophyll 52.4%, Epidermal 12.8%, Stele 31.2%, Meristematic 0%, Guard 0%
- Leiden clustering (res=0.5): 6 clusters, ARI=0.033, NMI=0.040 vs. CellType

### S4.6 Limitations
- `seq_length` was reduced from 1500 to 500 for CPU feasibility, which may degrade embedding quality.
- scPlantLLM was trained predominantly on root tissue; the atlas is whole-seedling, creating a domain shift.
- Zero-shot performance is expected to be lower than fine-tuned approaches; the results should be interpreted as a baseline for plant foundation model evaluation.

### S4.7 Code
```python
# scPlantLLM zero-shot inference (simplified)
import torch
from model import TransformerModel  # scPlantLLM architecture

# Model with n_input_bins=103 (NOT 101)
model = TransformerModel(ntoken=185622, d_model=512, nhead=8,
                         d_hid=512, nlayers=6, nlayers_cls=3,
                         n_cls=44, n_input_bins=103)

# Load converted weights
state_dict = torch.load('model_params/scPlantLLM_model_converted.pth', map_location='cpu')
model.load_state_dict(state_dict, strict=True)
model.eval()

# Inference for each cell
with torch.no_grad():
    for batch in dataloader:
        # In category mode: replace padding with N_BINS
        batch_val[batch_val == PAD_VALUE] = N_BINS
        output = model(batch_genes, batch_val, src_key_padding_mask=mask)
        pred = output.argmax(dim=1)
```

---

## S5. Clustering and Differential Expression

### S5.1 Clustering
scPlantLLM embeddings (512-dimensional) were extracted for all 9,998 cells. Leiden clustering was performed using scanpy with resolution 0.5, producing 6 clusters. Cluster quality was assessed by Adjusted Rand Index (ARI=0.033) and Normalized Mutual Information (NMI=0.040) against the atlas ground-truth cell type annotations.

### S5.2 Differential Expression
Differentially expressed genes (DEGs) were identified for each Leiden cluster using scanpy's `rank_genes_groups` function with the Wilcoxon rank-sum test. The top 20 DEGs per cluster were retained. Results were saved as cluster-specific gene markers.

### S5.3 Code
```python
import scanpy as sc
sc.pp.neighbors(adata, use_rep='X_scplantllm', n_neighbors=15)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
```

---

## S6. PlantCellChat Cell-Cell Communication Analysis

### S6.1 Database
PlantCellChat (https://github.com/mrliuw/PlantCellChat) uses a curated *Arabidopsis* ligand-receptor interaction database containing 3,140 interactions across 14 signaling pathways (Ca2+, ABA, BR, chitin, flg22, CTK, ET, GA, IAA, K+, Mg2+, NO3−, SA, elf18).

### S6.2 Analysis Parameters
- **Thresholds**: P-value threshold = 0.1, probability threshold = 0, fold-change threshold = 0 (relaxed from defaults to maximize discovery given the small cell sample)
- **Cell types**: 5 cell types (Epidermal, Mesophyll, Stele, Meristematic, Guard)
- **Communication strength**: Calculated as the sum of interaction probabilities across all ligand-receptor pairs for each source-target cell type pair.

### S6.3 Key Results
- Ca2+ signaling dominated with total strength 249.5, 4.6× the next pathway (BR=54.5)
- 371 unique Ca2+ ligand-receptor interactions measured
- Top Ca2+ LR pair: AT1G20620_AT5G19450 (strength=10.02)
- Communication strength matrix showed Meristematic→Meristematic as the strongest self-communication (92.0)

### S6.4 Limitations
- Guard cells had only 21 cells in the subsample, making CCC estimates for this cell type unstable.
- The atlas was grown under terrestrial conditions; CCC patterns may differ under spaceflight.
- Thresholds were relaxed from PlantCellChat defaults, which may increase false positive rates.

---

## S7. GO Enrichment Analysis

### S7.1 Method
Gene Ontology (GO) and KEGG pathway enrichment analysis was performed using the g:Profiler web service (Raudvere et al., Nucleic Acids Research, 2019) via the `gprofiler2` R package (organism: `athaliana`).

### S7.2 Parameters
- **Background**: Custom background of 2,000 genes (the top variable features from which LASSO selected), loaded from `all_feature_stability.csv`. This is the statistically correct background because LASSO selected from these 2,000 features, not the whole genome.
- **Correction**: Benjamini-Hochberg FDR (padj < 0.05)
- **Sources**: GO:BP (Biological Process), GO:MF (Molecular Function), GO:CC (Cellular Component), KEGG
- **Queries**: Three separate enrichment analyses were run:
  1. All 85 LASSO panel genes (combined)
  2. 41 flight-up genes (mean coefficient > 0)
  3. 44 ground-up genes (mean coefficient < 0)

### S7.3 Key Results
- **All 85 genes**: 2 significant terms — KEGG Protein processing in ER (p=2.6×10⁻⁸, 9 genes), GO:BP response to heat (p=0.011, 8 genes)
- **Flight-up (41 genes)**: 11 significant terms — response to heat (p=4.3×10⁻³), cellular response to hypoxia/oxygen levels, protein modification by small protein conjugation, KEGG protein processing in ER (p=1.3×10⁻¹¹, 9 genes)
- **Ground-up (44 genes)**: 1 significant term — flavonoid biosynthetic process (p=0.041, 4 genes)

The directed analysis reveals a clear biological split: spaceflight-upregulated genes are enriched for stress/heat-shock/ER-proteostasis pathways, while ground-upregulated genes are enriched for flavonoid/secondary metabolism.

### S7.4 Code
```r
library(gprofiler2)
background <- read.csv("all_feature_stability.csv")$feature  # 2000 genes
genes_all <- biomarker_panel$feature  # 85 genes
genes_pos <- biomarker_panel$feature[biomarker_panel$mean_coefficient > 0]  # 41
genes_neg <- biomarker_panel$feature[biomarker_panel$mean_coefficient < 0]  # 44

res <- gost(query = genes_all, organism = "athaliana",
            sources = c("GO:BP","GO:MF","GO:CC","KEGG"),
            significant = FALSE, correction_method = "fdr",
            custom_bg = background, evcodes = TRUE)
sig <- res$result[res$result$p_value < 0.05, ]
```

---

## S8. ggPlantmap Spatial Visualization

### S8.1 Method
Spatial expression visualization was performed using the ggPlantmap R package (Jo & Kajala, Journal of Experimental Botany, 2024), which provides anatomical template maps of *Arabidopsis* organs. Cell type-specific expression values were mapped onto tissue regions of interest (ROIs) using the package's merge and heatmap functions.

### S8.2 Cell Type to ROI Mapping

**Leaf cross-section:**
| Cell type | ROI |
|-----------|-----|
| Epidermal | epidermis.adaxial, epidermis.abaxial |
| Mesophyll | Parenchima.palisade, Parenchima.sponge |
| Stele | vascularbundle.xylem, vascularbundle.phloem, vascularbundle.bundlesheet |
| Guard | epidermis.stomata |

**Root tip cross-section:**
| Cell type | ROI |
|-----------|-----|
| Epidermal | Epidermis |
| Mesophyll | Cortex |
| Stele | Procambium, Xylem, Phloem, Pericycle |
| Meristematic | Columella |

### S8.3 Code
```r
library(ggPlantmap)
library(cowplot)
# Merge expression data with map
merged <- ggPlantmap.merge(map_data, values_df, "ROI.name")
# Heatmap with guide default gradient
p <- ggPlantmap.heatmap(merged, value.quant = value) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_void() + coord_fixed()
# Multi-panel composite
cowplot::plot_grid(p1, p2, ncol = 1, labels = c("a", "b"))
```

---

## S9. Ca2+ Pathway Network Visualization

### S9.1 Method
The Ca2+ signaling pathway was visualized as a directed network using igraph and ggraph. The network was constructed with 20 nodes representing key Ca2+ signaling components and 43 edges representing measured and literature-based interactions.

### S9.2 Node Selection
Nodes were selected based on two criteria:
1. **Data-driven**: Presence in the PlantCellChat Ca2+ ligand-receptor database with measured communication strength
2. **Literature-driven**: Known biological relevance to Ca2+ signaling and spaceflight stress response

Node types: Ca2+ ligand (1), CBL sensors (3: CBL1, CBL9, CBL2), CIPK kinases (3: CIPK23, CIPK1, CIPK9), CaM/CML sensors (5: CML38, CML42, CML24, CaM7, CaM3), CDPK kinases (2: CDPK3, CDPK4), downstream targets (4: CAX3, NRT2.1, SOS1, ANNAT1), LASSO biomarkers (2: HSP22.0, TRAF-like).

### S9.3 Edge Types
- **Primary (22 edges)**: Ca2+ → receptor interactions with measured CCC strength (solid blue, width proportional to strength)
- **Cascade (13 edges)**: CBL → CIPK → target phosphorylation/activation cascades from literature (dashed grey)
- **Branch (5 edges)**: CaM/CDPK → downstream signaling branches (dotted pink)
- **LASSO links (3 edges)**: Connections between Ca2+ signaling and LASSO biomarker genes (dotted orange)

### S9.4 Layout
Node placement was computed using the stress-majorization algorithm (`graphlayouts::layout_with_stress()`), which minimizes the P-stress across all node pairs. The Ca2+ ligand node was manually shifted to the top of the layout for visual clarity. Node size is proportional to mean expression in the atlas; node color indicates functional type.

### S9.5 Code
```r
library(igraph)
library(ggraph)
library(graphlayouts)
# Build network
ca2_network <- graph_from_data_frame(d = edge_table, vertices = node_table, directed = TRUE)
# Stress-majorization layout
set.seed(42)
layout_stress <- layout_with_stress(ca2_network)
# Plot with ggraph
p <- ggraph(ca2_network, layout = "manual", x = layout_stress[,1], y = layout_stress[,2]) +
  geom_edge_link(aes(filter = edge_type == "primary", width = strength), ...) +
  geom_node_point(aes(fill = type, size = expression), shape = 21) +
  geom_node_text(aes(label = alias), repel = TRUE) +
  theme_void()
```

---

## S10. Gene Function Annotation

### S10.1 Method
Gene function annotations for the top 4 LASSO biomarker genes were compiled from multiple databases:
- **UniProt**: Protein names, families, and molecular weights
- **TAIR/Arabidopsis.org**: Gene aliases and functional descriptions
- **SUBAcon**: Subcellular localization predictions
- **KEGG**: Pathway assignments
- **Gene Ontology**: Molecular function, biological process, cellular component
- **Literature**: Spaceflight relevance and known stress response functions

### S10.2 Key Findings
- **AT4G10250 (HSP22.0)**: HSP20/alpha-crystallin family chaperone, ER and chloroplast localized, ABA-auxin crosstalk, flight-enriched positive predictor (coef=+0.723)
- **AT3G07365**: Natural antisense transcript (NAT), non-coding RNA overlapping AT3G46230, strongest positive predictor (coef=+0.828), highlights non-coding RNA regulation in spaceflight
- **AT2G14247**: Small chloroplast-localized protein (78 aa, 8.65 kDa), function unknown, ground-enriched negative predictor (coef=−0.364)
- **AT4G01390**: TRAF-like/MATH domain protein, signal transduction, senescence-associated, H2O2-responsive, strongest negative predictor (coef=−0.819)

---

## Software Versions

| Software | Version |
|----------|---------|
| Python | 3.11+ |
| scanpy | 1.11.5 |
| pandas | 2.3.3 |
| matplotlib | 3.11.0 |
| scikit-learn | (for LASSO) |
| PyTorch | (for scPlantLLM) |
| R | 4.4.x |
| Seurat | (for atlas loading) |
| glmnet | (for LASSO) |
| igraph | (for network) |
| ggraph | (for network visualization) |
| graphlayouts | (for stress-majorization layout) |
| ggPlantmap | (from GitHub, Jo & Kajala 2024) |
| gprofiler2 | (for GO enrichment) |
| PlantCellChat | (from GitHub, mrliuw) |
| ggkegg | (from GitHub, noriakis; Figure 8 panel a) |
| KEGGREST | (Bioconductor; KEGG pathway retrieval) |
| png | (for KEGG raster I/O) |
| plotly | 6.8.0 (for Sankey diagram) |
| kaleido | 0.2.1 (for static image export) |
| python-docx | (for manuscript generation) |
| pandoc | (for PDF conversion) |
| weasyprint | (for PDF rendering) |

## Data Availability

- **OSDR data**: Available from https://visualization.osdr.nasa.gov/biodata/api/ (study IDs: OSD-37, OSD-678, OSD-38, OSD-321, OSD-120, OSD-624)
- **Seedling atlas**: GEO accession GSE226097 (Lee et al., Nature Plants 2025)
- **scPlantLLM pretrained weights**: Available from the scPlantLLM repository (Cao et al., 2025)
- **All analysis code**: Available in the companion repository (see Data Availability in main manuscript)

## S11. Ca2+/K+ Crosstalk Circuit Analysis

### S11.1 Network Expansion
The Ca2+ pathway network (S9) was expanded to include six major K+ channels identified from the PlantCellChat K+ pathway database:
- **AKT1** (AT2G26650) — primary Ca2+-regulated K+ uptake channel
- **AKT2** (AT5G46240) — phloem K+ channel
- **GORK** (AT4G18290) — guard cell outward-rectifying K+ channel
- **KAT1** (AT5G37500) — guard cell inward-rectifying K+ channel
- **KAT2** (AT4G22200) — K+ channel
- **KC1** (AT4G30960) — K+ channel subunit, modulates AKT1

The CBL9-CIPK23-AKT1 crosstalk edge was added based on literature evidence [12,13,15]. CBL9 (AT5G24270) is present in both the Ca2+ and K+ PlantCellChat pathway databases, serving as the molecular bridge between the two signaling systems. The resulting network contains 27 nodes and 67 edges (28 Ca2+ primary, 12 K+ primary, 14 cascade, 5 branch, 4 K+ channel, 1 crosstalk, 3 LASSO).

### S11.2 Co-expression Analysis
Pearson correlation was computed for all 25 Ca2+/K+ circuit genes across:
- All cells (global, 41,314 cells)
- Per cell type (Epidermal, Mesophyll, Stele, Meristematic, Guard)

Key finding: Guard cells (n=89) show the strongest co-expression module:
- CDPK4 ~ ANNAT1: r = 0.559
- CIPK23 ~ CaM3: r = 0.460
- CDPK3 ~ GORK: r = 0.382
- CBL9 ~ AKT1: r = 0.337

All other cell types show weak correlations (r < 0.1), suggesting the circuit operates primarily in guard cells.

### S11.3 Cell-Cell Communication Circuit
Ca2+ and K+ communication strengths were extracted from PlantCellChat outputs and aggregated into 5×5 source-target matrices. Key findings:
- Ca2+ total strength: 249.5 (10.2-fold higher than K+ at 24.5)
- Guard cells: dominant Ca2+ source (outgoing total 57.1; top target: Meristematic 12.14)
- Epidermal cells: primary K+ source (outgoing total 6.06)
- Top K+ LR pair: KC1 (AT4G30960) → AKT1 (strength 5.01)

### S11.4 Organ Marker Inference
Marker-based organ inference was attempted using 15 tissue-specific markers:
- Root: SCR, PLT1, PLT2, WOL, CPC
- Leaf: LHCB1.1, FSD1, RBCS, COR15A
- Hypocotyl: ATHB-1, ATHB-7
- Seed: LEC2, 2S3
- Guard: KAT1_gc, MYB60

Result: Marker-based inference was NOT reliable. Photosynthetic markers (RBCS 93-98%, FSD1 72-85%) were ubiquitous across all cell types. Root markers (SCR <3%, PLT1 <11%, WOL <1.2%) were near-absent. Seed markers (LEC2 <0.2%, 2S3 <0.6%) were near-zero (expected for 6-day seedling). Cell-type proxy mapping was used instead.

### S11.5 Figure 8 Construction (KEGG overlay + molecular cascade)
Figure 8 is a two-panel figure summarizing the CBL9-CIPK23-AKT1 cascade:
- **Panel A (KEGG pathway overlay)**: The cascade was mapped onto KEGG reference pathway `ath04075` (Plant hormone signal transduction) using the `ggkegg` R package (Sato et al., GitHub: noriakis/ggkegg). The KEGG pathway image was retrieved via `KEGGREST::keggGet('ath04075', 'image')` and rendered as a background raster with `ggplot2::annotation_raster()`. Circuit genes present in KEGG were highlighted on their native nodes: CDPK4 (AT4G09570), CDPK3 (AT4G23650), and CIPK23 (AT4G35310) share one KEGG node (KO:K13412; orange), and CaM3 (AT3G56800) and CaM7 (AT3G43810) share another (KO:K02183; pink). Circuit genes absent from any KEGG Arabidopsis pathway (CBL9, CBL1, CBL2, CML24, AKT1, KC1, ANNAT1) were added as custom nodes colored by guard cell mean expression (white→blue gradient). Cascade arrows trace Ca2+ → CBL9 → CIPK23 → AKT1 → K+ uptake, with guard cell co-expression r values labeled on key edges. Rendered with `ggplot2` at 300 dpi.
- **Panel B (molecular cascade)**: A Plotly Sankey diagram. Ca2+ → CBL1/CBL9/CaM3/CML24 (sensors) → CIPK23/CDPK4/CDPK3 (kinases) → AKT1/ANNAT1/KC1 (effectors) → K+ uptake. Flow width = guard cell mean expression × 100. Edge labels show guard cell co-expression r values. Static export via kaleido (v0.2.1) at scale=3.

The two panels were composited side-by-side with panel labels using PIL (Python). Rationale for the KEGG overlay: it places the cascade in the context of the canonical plant hormone signaling map, making explicit which components are curated KEGG entities (CDPK/CIPK kinases, calmodulins) versus circuit-specific genes (CBL sensors, K+ channels) that KEGG does not resolve for Arabidopsis.

### S11.6 ggPlantmap Montage
Three ggPlantmap composites were generated and combined into a supplementary montage:
1. **Seedling map** (ggPm.At.seedling.saltdrought): cotyledon, hypocotyl, root regions
2. **Root tip longitudinal** (ggPm.At.roottip.longitudinal): 12 root zones
3. **Leaf cross-section** (ggPm.At.leaf.crosssection): 8 leaf regions

Cell-type expression was mapped to organ regions as proxy:
- Cotyledon = Mesophyll + Guard average
- Hypocotyl = Epidermal + Stele average
- Root = Stele + Meristematic average

### S11.7 Code
```r
# Panel A: ggKEGG cascade overlay on ath04075 (R)
library(ggkegg); library(KEGGREST); library(ggplot2); library(png)

kegg_img <- keggGet("ath04075", "image")   # 3D raster array [h, w, 4]
img_w <- ncol(kegg_img); img_h <- nrow(kegg_img)

# Custom nodes for genes not in KEGG, colored by guard cell expression
p <- ggplot() +
  annotation_raster(kegg_img, xmin=0, xmax=img_w, ymin=0, ymax=img_h) +
  # cascade arrows Ca2+ -> CBL9 -> CIPK23 -> AKT1 -> K+ uptake
  geom_segment(data=cascade, aes(x=x_from, y=y_from, xend=x_to, yend=y_to),
               arrow=arrow(type="closed")) +
  # highlight native KEGG nodes (CDPK/CIPK, CaM) + custom overlay nodes
  geom_point(data=non_kegg, aes(x=plot_x, y=plot_y, fill=guard_expr),
             shape=21, size=8) +
  scale_fill_gradient(low="white", high="#0279EE") +
  coord_cartesian(xlim=c(0,img_w), ylim=c(0,img_h)) + scale_y_reverse() +
  theme_void()
ggsave("figure8a_ggkegg.png", p, width=10, height=15, dpi=300)
```

```python
# Panel B: molecular cascade Sankey (Plotly) + composite (PIL)
import plotly.graph_objects as go
from PIL import Image

sankey_b = go.Sankey(
    node=dict(label=cascade_nodes, color=node_colors_b),
    link=dict(source=cascade_src, target=cascade_tgt,
              value=cascade_val, color=cascade_color, label=r_labels))
fig_b = go.Figure(data=[sankey_b])
fig_b.write_image('figure8b_sankey.png', scale=3)

# Composite panel A (KEGG) + panel B (Sankey) side-by-side
a = Image.open('figure8a_ggkegg.png'); b = Image.open('figure8b_sankey.png')
# ... resize to matching height, paste side-by-side, add panel labels ...
composite.save('Figure8_sankey_ca2_k_cascade.png')
```

```r
# ggPlantmap montage (R + Python compositing)
# Individual maps generated in R (see script 12_ggplantmap_ca2_k_circuit.R)
# Montage assembled in Python using matplotlib
```

---

## S12. GO Enrichment Analysis (Detailed)

### S12.1 Rationale
The 85-gene LASSO panel was selected from 2,000 variable features. The statistically correct background for enrichment is therefore these 2,000 features (not the whole genome), because LASSO only had access to these features during selection.

### S12.2 Three Directed Queries
The panel was split by coefficient direction to identify biological processes specific to flight-upregulated versus ground-upregulated genes:
1. **All 85 genes** (combined panel)
2. **41 flight-up genes** (mean coefficient > 0) — genes whose expression increases in spaceflight
3. **44 ground-up genes** (mean coefficient < 0) — genes whose expression decreases in spaceflight

### S12.3 Results Summary

| Query | Significant terms | Top term | p-value |
|-------|------------------|----------|---------|
| All 85 | 2 | KEGG Protein processing in ER | 2.6×10⁻⁸ |
| Flight-up (41) | 11 | KEGG Protein processing in ER | 1.3×10⁻¹¹ |
| Ground-up (44) | 1 | Flavonoid biosynthetic process | 0.041 |

### S12.4 Biological Interpretation
The directed analysis reveals a clear biological split:
- **Flight-up genes**: enriched for stress/heat-shock/ER-proteostasis pathways (response to heat, cellular response to hypoxia, protein modification by small protein conjugation, protein processing in ER)
- **Ground-up genes**: enriched for flavonoid/secondary metabolism (flavonoid biosynthetic process)

This split is consistent with the known spaceflight response: upregulation of stress machinery (heat shock proteins, ER folding capacity) and downregulation of secondary metabolism (energy reprioritization).

### S12.5 Visualization
Two figures were generated:
- **Dot plot** (Supplementary Figure S1): significant terms across all 3 queries, dot size = gene count, color = -log10(p)
- **Directed bar plot** (Supplementary Figure S2): -log10(p) separated by direction (flight-up vs ground-up)

### S12.6 Code
```r
library(gprofiler2)
background <- read.csv("all_feature_stability.csv")$feature  # 2000 genes
genes_all <- biomarker_panel$feature  # 85 genes
genes_pos <- biomarker_panel$feature[biomarker_panel$mean_coefficient > 0]  # 41
genes_neg <- biomarker_panel$feature[biomarker_panel$mean_coefficient < 0]  # 44

# Run 3 queries
for (query_name in c("All_85", "Flight_up_41", "Ground_up_44")) {
  if (query_name == "All_85") query_genes <- genes_all
  else if (query_name == "Flight_up_41") query_genes <- genes_pos
  else query_genes <- genes_neg

  res <- gost(query = query_genes, organism = "athaliana",
              sources = c("GO:BP","GO:MF","GO:CC","KEGG"),
              significant = FALSE, correction_method = "fdr",
              custom_bg = background, evcodes = TRUE)
  sig <- res$result[res$result$p_value < 0.05, ]
  write.csv(sig, paste0("go_enrichment_", query_name, ".csv"))
}
```
