# Quality control for NanoString GeoMx DSP spatial transcriptomics data
# Framework: GeoMxTools
# QC performed at the ROI (AOI) level

library(GeoMxTools)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)

# --------------------------------------------------
# Load imported GeoMx object
# --------------------------------------------------
geomx_obj <- readRDS("geomx_imported_raw.rds")

# --------------------------------------------------
# Extract ROI-level QC metrics
# --------------------------------------------------
roi_qc <- pData(geomx_obj) %>%
  as_tibble() %>%
  select(
    ROI = roi_id,
    Segment = segment,
    Area = area,
    Nuclei = nuclei,
    LibrarySize = total_counts
  )

# --------------------------------------------------
# Visual QC assessment
# --------------------------------------------------
ggplot(roi_qc, aes(x = LibrarySize)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Library size distribution across ROIs")

ggplot(roi_qc, aes(x = Nuclei, y = LibrarySize)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Library size vs nuclei count",
    x = "Nuclei count",
    y = "Library size"
  )

# --------------------------------------------------
# Identify low-quality ROIs
# --------------------------------------------------
low_lib_threshold <- quantile(roi_qc$LibrarySize, 0.05, na.rm = TRUE)

roi_qc <- roi_qc %>%
  mutate(
    LowQuality = LibrarySize < low_lib_threshold
  )

# --------------------------------------------------
# Filter ROIs conservatively
# --------------------------------------------------
high_quality_rois <- roi_qc %>%
  filter(!LowQuality) %>%
  pull(ROI)

geomx_qc <- geomx_obj[, pData(geomx_obj)$roi_id %in% high_quality_rois]

# --------------------------------------------------
# Save QC-filtered object
# --------------------------------------------------
saveRDS(geomx_qc, "geomx_qc_filtered.rds")
