# Normalization and feature selection for spatial transcriptomics

library(Seurat)

spatial_obj <- readRDS("spatial_qc_filtered.rds")

# -----------------------------
# Log-normalization
# -----------------------------
spatial_obj <- NormalizeData(
  spatial_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# -----------------------------
# Identify variable features
# -----------------------------
spatial_obj <- FindVariableFeatures(
  spatial_obj,
  selection.method = "vst",
  nfeatures = 2000
)

saveRDS(spatial_obj, "spatial_normalized.rds")
