\# Methods Appendix: Spatial Omics Workflows



\## Overview



This document describes a complete spatial transcriptomics analysis

workflow, from raw data import to spatial visualization, differential

expression, and pathway-level interpretation. The workflows are designed

for NanoString GeoMx DSP and similar spatial omics platforms.



The methods emphasize reproducibility, analytical rationale, and

clarity, rather than serving as a turnkey software package.



---



\## 0. Software Environment and Setup



Analyses are performed using R (≥4.1) and Bioconductor packages suited

for spatial transcriptomics. Key dependencies include Seurat,

GeoMxTools, SpatialExperiment, and tidyverse-based visualization tools.

All software versions and reference files are documented to support

reproducibility.



---



\## 1. Data Import and Preprocessing



Spatial transcriptomics data are imported from platform-specific output

formats, including raw count matrices and region-level metadata. Spatial

coordinates, imaging metadata, and experimental annotations are merged

with expression data to create a unified analysis object.



Raw counts are preserved to allow transparent downstream processing and

re-analysis.



---



\## 2. Quality Control of Spatial Regions



Quality control is performed at the region-of-interest (ROI) level.

Metrics include total counts, number of detected genes, and platform-

specific control features. Outlier regions are identified using

visualization and summary statistics and may be excluded if they fail

quality thresholds.



QC decisions are documented to ensure interpretability and reproducibility.



---



\## 3. Normalization and Feature Selection



Expression data are normalized using methods appropriate for spatial

omics platforms, including log-normalization or model-based approaches.

Highly variable genes are identified to support dimensionality

reduction and clustering.



Normalization choices are guided by platform characteristics and data

distribution.



---



\## 4. Dimensionality Reduction and Integration



Dimensionality reduction methods such as PCA and UMAP are applied to

normalized data to explore global structure and spatial heterogeneity.

Where multiple samples or batches are present, integration strategies

are used to mitigate technical variation.



---



\## 5. Clustering and Spatial Domain Annotation



Spatial regions are clustered using graph-based approaches. Clusters

are annotated using canonical markers, spatial context, and known tissue

structure. Annotations are assigned conservatively and may be updated as

additional biological validation becomes available.



---



\## 6. Differential Expression Analysis



Differential expression testing is performed between spatial domains or

conditions using appropriate statistical models. Covariates such as

batch or region-level effects may be included where applicable.

Differential expression results are used for downstream biological

interpretation.



---



\## 7. Spatial Visualization and Mapping



Spatial expression patterns are visualized using coordinate-based plots

that map gene expression and cluster identity back onto tissue context.

These visualizations support interpretation of spatial organization and

biological structure.



---



\## 8. Pathway Enrichment Analysis



Pathway-level interpretation is performed using both gene set

enrichment analysis (GSEA) and over-representation analysis (ORA).

Curated gene sets are used to identify biological processes associated

with spatial domains or conditions.



---



\## 9. Reporting and PI-facing Outputs



Results are summarized in concise tables and figures designed for rapid

review by principal investigators. Outputs include ranked pathway

tables, spatial maps, and annotated figures, exported in spreadsheet-

compatible formats for downstream use.



