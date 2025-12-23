# Spatial visualization of gene expression and clusters

library(Seurat)
library(patchwork)

spatial_obj <- readRDS("spatial_clustered.rds")

# -----------------------------
# UMAP cluster visualization
# -----------------------------
DimPlot(spatial_obj, reduction = "umap", label = TRUE)

# -----------------------------
# Feature plots (examples)
# -----------------------------
FeaturePlot(
  spatial_obj,
  features = head(VariableFeatures(spatial_obj), 4),
  ncol = 2
)
