# Clustering spatial domains

library(Seurat)

spatial_obj <- readRDS("spatial_umap.rds")

# -----------------------------
# Graph-based clustering
# -----------------------------
spatial_obj <- FindNeighbors(spatial_obj, dims = 1:20)
spatial_obj <- FindClusters(spatial_obj, resolution = 0.4)

# -----------------------------
# Visualization
# -----------------------------
DimPlot(spatial_obj, reduction = "umap", label = TRUE)

saveRDS(spatial_obj, "spatial_clustered.rds")
