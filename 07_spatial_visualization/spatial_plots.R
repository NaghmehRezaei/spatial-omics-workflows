# Spatial visualization for NanoString GeoMx DSP data
# Purpose: map expression and annotations back to tissue context

library(GeoMxTools)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)

# --------------------------------------------------
# Load QC-filtered and normalized GeoMx object
# --------------------------------------------------
geomx_norm <- readRDS("geomx_normalized.rds")

# --------------------------------------------------
# Extract spatial metadata
# --------------------------------------------------
meta <- pData(geomx_norm) %>%
  as_tibble()

# Expected metadata columns (example):
# x, y         -> spatial coordinates
# Group        -> spatial domain / condition
# roi_id       -> ROI identifier

# --------------------------------------------------
# Plot spatial distribution of ROIs
# --------------------------------------------------
ggplot(meta, aes(x = x, y = y, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Spatial distribution of regions of interest",
    x = "X coordinate",
    y = "Y coordinate"
  )

# --------------------------------------------------
# Visualize expression of selected genes
# --------------------------------------------------
expr <- assay(geomx_norm, "logCPM")

genes_to_plot <- rownames(expr)[1:3]  # example subset

expr_long <- expr[genes_to_plot, ] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(meta) %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "Gene",
    values_to = "Expression"
  )

ggplot(expr_long, aes(x = x, y = y, color = Expression)) +
  geom_point(size = 3) +
  facet_wrap(~ Gene) +
  coord_fixed() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Spatial gene expression patterns")
