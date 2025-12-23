# Import spatial transcriptomics data (NanoString GeoMx-style)
# Input: raw count matrix and ROI-level metadata
# Output: unified Seurat object for spatial analysis

library(Seurat)
library(tidyverse)

# -----------------------------
# Load expression data
# -----------------------------
# Rows: genes
# Columns: regions of interest (ROIs)
counts <- read_tsv("data/spatial_counts.tsv")

counts_mat <- counts %>%
  column_to_rownames(var = colnames(counts)[1]) %>%
  as.matrix()

# -----------------------------
# Load ROI metadata
# -----------------------------
roi_meta <- read_tsv("data/roi_metadata.tsv")

# -----------------------------
# Create Seurat object
# -----------------------------
spatial_obj <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = roi_meta
)

saveRDS(spatial_obj, "spatial_raw_seurat.rds")
