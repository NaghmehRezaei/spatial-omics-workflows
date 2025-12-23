# 🧬 Spatial Omics Workflows

This repository documents end-to-end **spatial transcriptomics analysis
workflows**, with a focus on **NanoString GeoMx DSP** and related spatial
omics platforms.

📄 **Complete analytical methods are documented in [`METHODS.md`](METHODS.md).**

---

## Overview

The workflows in this repository reflect practical, analysis-driven
approaches used for spatial transcriptomics data, with an emphasis on
region-level modeling, biological interpretability, and reproducibility.
They are written in a **methods-appendix style**, mirroring how spatial
omics analyses are typically described in publications and grant
supplements.

---

## Workflow summary

GeoMx DSP counts and ROI metadata  
→ ROI-level quality control (library size, nuclei, control probes)  
→ Q3 / logCPM normalization  
→ limma-voom differential expression analysis  
→ exploratory PCA and UMAP (visualization only)  
→ spatial expression mapping across tissue context  
→ gene set enrichment analysis (GSEA)  
→ PI-facing summary tables and figures

---

## Analytical framework

Spatial transcriptomics analyses are performed using a combination of
**GeoMxTools** for NanoString DSP–specific data handling and quality
control, **limma/voom** for region-level differential expression
analysis, and **fgsea** for pathway enrichment. **Seurat** is used for
dimensionality reduction and visualization, but not for statistical
modeling.

This framework reflects common analytical practice across GeoMx DSP
studies rather than a single rigid pipeline.

---

## Scope of analysis

The workflows cover:

- Import of NanoString GeoMx DSP expression data and ROI metadata  
- Quality control at the region-of-interest (ROI) level  
- Normalization and feature filtering appropriate for spatial data  
- Dimensionality reduction and exploratory visualization  
- Identification and annotation of spatial domains  
- Differential expression between spatial regions or conditions  
- Pathway enrichment analysis (GSEA and ORA)  
- Generation of PI-facing summary tables and figures  

---

## Intended audience

- Principal investigators and senior scientists  
- Computational biology collaborators  
- Bioinformatics reviewers  
- Trainees seeking structured spatial omics workflows  

---

## Design philosophy

This repository prioritizes transparent, statistically sound analysis
over black-box automation. Each analytical step reflects common practice
in NanoString GeoMx DSP studies and is designed to be inspected, adapted,
and interpreted rather than executed as a fixed pipeline.

---

## Reproducibility and data availability

Raw spatial data and sensitive metadata are intentionally excluded.
Scripts are designed to be adapted to specific experimental designs,
platform configurations, and institutional computing environments.
