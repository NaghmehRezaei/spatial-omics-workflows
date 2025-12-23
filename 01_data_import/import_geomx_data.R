# Import and preprocess NanoString GeoMx DSP spatial transcriptomics data
# Platform: NanoString GeoMx DSP
# Framework: GeoMxTools
# Output: GeoMxTools-compatible object for downstream analysis

library(GeoMxTools)
library(SummarizedExperiment)
library(tidyverse)

# --------------------------------------------------
# Load GeoMx DSP data
# --------------------------------------------------
# Expected input: GeoMxTools-compatible object containing
# expression counts, ROI metadata, and control probes

geomx_obj <- readRDS("data/geomx_raw_data.rds")

# --------------------------------------------------
# Basic object checks
# --------------------------------------------------
stopifnot(
  inherits(geomx_obj, "NanoStringGeoMxSet") ||
    inherits(geomx_obj, "SummarizedExperiment")
)

# --------------------------------------------------
# Preserve raw counts
# --------------------------------------------------
assayNames(geomx_obj)

# Raw counts are retained to allow transparent
# normalization and reprocessing downstream

saveRDS(geomx_obj, "geomx_imported_raw.rds")
