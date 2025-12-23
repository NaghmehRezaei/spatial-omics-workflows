# Quality control for spatial transcriptomics ROIs
# Output: QC-filtered Seurat object

library(Seurat)
library(ggplot2)

spatial_obj <- readRDS("spatial_raw_seurat.rds")

# -----------------------------
# QC metrics
# -----------------------------
spatial_obj[["nFeature_RNA"]] <- spatial_obj$nFeature_RNA
spatial_obj[["nCount_RNA"]] <- spatial_obj$nCount_RNA

# -----------------------------
# Visual QC
# -----------------------------
VlnPlot(
  spatial_obj,
  features = c("nFeature_RNA", "nCount_RNA"),
  pt.size = 0.1
)

# -----------------------------
# Conservative filtering
# -----------------------------
spatial_obj <- subset(
  spatial_obj,
  subset = nFeature_RNA > 200 & nCount_RNA > 500
)

saveRDS(spatial_obj, "spatial_qc_filtered.rds")
