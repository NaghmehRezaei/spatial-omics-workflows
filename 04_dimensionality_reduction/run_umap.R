# Dimensionality reduction for NanoString GeoMx DSP data
# Purpose: exploratory visualization (not inferential analysis)
# Methods: PCA followed by UMAP

library(GeoMxTools)
library(SummarizedExperiment)
library(Seurat)
library(tidyverse)

# --------------------------------------------------
# Load normalized GeoMx object
# --------------------------------------------------
geomx_norm <- readRDS("geomx_normalized.rds")

# --------------------------------------------------
# Extract normalized expression (logCPM)
# --------------------------------------------------
expr <- assay(geomx_norm, "logCPM")

# --------------------------------------------------
# Create Seurat object for visualization only
# --------------------------------------------------
seurat_obj <- CreateSeuratObject(
  counts = expr,
  meta.data = pData(geomx_norm)
)

# --------------------------------------------------
# PCA
# --------------------------------------------------
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# --------------------------------------------------
# UMAP
# --------------------------------------------------
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# --------------------------------------------------
# Visualization
# --------------------------------------------------
DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "Group",
  label = TRUE
)

saveRDS(seurat_obj, "geomx_umap_seurat.rds")
