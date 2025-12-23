# Normalization of NanoString GeoMx DSP spatial transcriptomics data
# Framework: GeoMxTools
# Strategy: Q3-style normalization followed by log-transformation

library(GeoMxTools)
library(SummarizedExperiment)
library(tidyverse)
library(edgeR)

# --------------------------------------------------
# Load QC-filtered GeoMx object
# --------------------------------------------------
geomx_qc <- readRDS("geomx_qc_filtered.rds")

# --------------------------------------------------
# Extract raw count matrix
# --------------------------------------------------
counts <- assay(geomx_qc, "counts")

# --------------------------------------------------
# Calculate library-size–based normalization factors
# --------------------------------------------------
dge <- DGEList(counts = counts)

dge <- calcNormFactors(
  dge,
  method = "upperquartile"  # Q3 normalization
)

# --------------------------------------------------
# Log-CPM transformation
# --------------------------------------------------
logcpm <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

# --------------------------------------------------
# Store normalized expression
# --------------------------------------------------
assay(geomx_qc, "logCPM") <- logcpm

saveRDS(geomx_qc, "geomx_normalized.rds")
