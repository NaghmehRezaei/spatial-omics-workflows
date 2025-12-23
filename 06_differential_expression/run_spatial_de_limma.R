# Differential expression analysis for NanoString GeoMx DSP data
# Method: limma with voom transformation
# Level: ROI-level modeling

library(GeoMxTools)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(tidyverse)

# --------------------------------------------------
# Load normalized GeoMx object
# --------------------------------------------------
geomx_norm <- readRDS("geomx_normalized.rds")

# --------------------------------------------------
# Extract expression matrix (logCPM)
# --------------------------------------------------
expr <- assay(geomx_norm, "logCPM")

# --------------------------------------------------
# Define experimental design
# --------------------------------------------------
# Example grouping variable: spatial domain or condition
meta <- pData(geomx_norm)

group <- factor(meta$Group)  # e.g., "Zone1", "Zone2", "Control", "Treatment"

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# --------------------------------------------------
# voom transformation
# --------------------------------------------------
v <- voom(expr, design, plot = TRUE)

# --------------------------------------------------
# Linear modeling
# --------------------------------------------------
fit <- lmFit(v, design)

# Example contrast (adjust as needed)
contrast_matrix <- makeContrasts(
  Contrast1 = levels(group)[2] - levels(group)[1],
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# --------------------------------------------------
# Extract DE results
# --------------------------------------------------
de_results <- topTable(
  fit2,
  coef = "Contrast1",
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Gene")

# --------------------------------------------------
# Save results
# --------------------------------------------------
write_tsv(de_results, "geomx_DE_limma_results.tsv")
