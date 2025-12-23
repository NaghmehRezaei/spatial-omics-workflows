# Methods Appendix: Spatial Omics Workflows

## Overview

This document describes a complete workflow for the analysis of spatial
transcriptomics data generated using **NanoString GeoMx DSP** and related
platforms. The methods are written to emphasize analytical rationale,
reproducibility, and biological interpretability, and are intended to
reflect practical workflows commonly used in spatial omics studies.

The pipeline is modular and adaptable, rather than a fixed turnkey
solution, allowing analytical choices to be tailored to specific
datasets and experimental designs.

---

## 0. Software environment and analytical framework

All analyses are performed in R using **GeoMxTools** for NanoString
GeoMx DSP–specific data handling and quality control. Region-level
differential expression analysis is conducted using **limma/voom** to
appropriately model variance in spatial transcriptomics data. Pathway
enrichment analysis is performed using **fgsea**.

**Seurat** is used exclusively for dimensionality reduction and
visualization and is not used for statistical modeling of differential
expression.

Specific normalization strategies, quality control thresholds, and
statistical parameters may be adapted based on platform characteristics
and experimental design.

---

## 1. Data import and preprocessing

Spatial transcriptomics data are imported from NanoString GeoMx DSP
outputs, including raw count matrices and region-of-interest (ROI)
metadata. Expression data, spatial coordinates, and experimental
annotations are merged to create a unified analysis object suitable for
downstream modeling.

Raw counts are preserved during import to allow transparent
reprocessing and downstream normalization.

---

## 2. Quality control of spatial regions

Quality control is performed at the ROI level using standard spatial
omics metrics, including library size, detected gene counts, and
platform-specific control probes. ROIs exhibiting outlier behavior are
identified through summary statistics and visualization and may be
excluded if they fail predefined quality thresholds.

All QC decisions are documented to ensure interpretability and
reproducibility.

---

## 3. Normalization and feature filtering

Expression data are normalized using methods appropriate for GeoMx DSP
data, such as Q3 normalization or log-transformed counts per million.
Feature filtering is applied to remove low-information genes prior to
downstream analyses.

Normalization choices are guided by data distribution and platform
characteristics rather than fixed assumptions.

---

## 4. Dimensionality reduction and exploratory analysis

Dimensionality reduction techniques, including principal component
analysis (PCA) and uniform manifold approximation and projection (UMAP),
are applied to normalized expression data to explore global structure
and spatial heterogeneity across ROIs.

These analyses are used for exploratory purposes and quality assessment
rather than inferential testing.

---

## 5. Spatial domain identification and annotation

ROIs are grouped into spatial domains using clustering approaches
informed by reduced-dimensional embeddings. Identified domains are
annotated based on known marker genes, spatial context, and tissue
architecture.

Annotations are assigned conservatively and may be updated as additional
biological validation becomes available.

---

## 6. Differential expression analysis

Differential expression analysis is performed between spatial domains or
experimental conditions using **limma/voom**. Where applicable,
covariates such as batch effects or sample-level variation are included
in the model design.

Differential expression results form the basis for downstream biological
interpretation and pathway analysis.

---

## 7. Spatial visualization and mapping

Gene expression patterns, spatial domains, and experimental conditions
are visualized using coordinate-based plots that map molecular signals
back onto tissue context. These visualizations support interpretation of
spatial organization and biological structure within the tissue.

---

## 8. Pathway enrichment analysis

Pathway-level interpretation is performed using both gene set enrichment
analysis (GSEA) and over-representation analysis (ORA). Curated gene
sets are used to identify biological processes associated with spatial
domains or experimental conditions, with emphasis on rank-based
approaches appropriate for spatial transcriptomics data.

---

## 9. Reporting and PI-facing outputs

Final results are summarized in concise tables and figures designed to
support rapid review by principal investigators. Outputs include ranked
pathway summary tables, spatial visualizations, and annotated figures,
exported in spreadsheet-compatible formats for downstream reporting,
presentation, and collaborative analysis.
