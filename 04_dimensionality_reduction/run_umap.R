# Dimensionality reduction for spatial transcriptomics

library(Seurat)

spatial_obj <- readRDS("spatial_normalized.rds")

# -----------------------------
# Scale and PCA
# -----------------------------
spatial_obj <- ScaleData(spatial_obj, features = VariableFeatures(spatial_obj))
spatial_obj <- RunPCA(spatial_obj, features = VariableFeatures(spatial_obj))

# -----------------------------
# UMAP embedding
# -----------------------------
spatial_obj <- RunUMAP(spatial_obj, dims = 1:20)

saveRDS(spatial_obj, "spatial_umap.rds")
