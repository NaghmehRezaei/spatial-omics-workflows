# Differential expression between spatial clusters

library(Seurat)
library(dplyr)

spatial_obj <- readRDS("spatial_clustered.rds")

Idents(spatial_obj) <- "seurat_clusters"

# -----------------------------
# Differential expression
# -----------------------------
markers <- FindAllMarkers(
  spatial_obj,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write_tsv(markers, "spatial_cluster_markers.tsv")
