# Install required packages for NanoString GeoMx spatial transcriptomics workflows

# --------------------------------------------------
# CRAN packages
# --------------------------------------------------
install.packages(c(
  "tidyverse",
  "ggplot2",
  "patchwork",
  "Seurat"
))

# --------------------------------------------------
# Bioconductor packages
# --------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GeoMxTools",
  "SummarizedExperiment",
  "SpatialExperiment",
  "limma",
  "edgeR",
  "fgsea",
  "msigdbr"
))
