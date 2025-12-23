# Pathway enrichment analysis for NanoString GeoMx DSP data
# Method: rank-based Gene Set Enrichment Analysis (GSEA)
# Input: limma differential expression results

library(readr)
library(dplyr)
library(fgsea)
library(msigdbr)

# --------------------------------------------------
# Load differential expression results
# --------------------------------------------------
de <- read_tsv("geomx_DE_limma_results.tsv")

# --------------------------------------------------
# Prepare ranked gene list
# --------------------------------------------------
# Ranking by moderated t-statistic is standard for limma-based GSEA
ranks <- de$t
names(ranks) <- de$Gene

ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

# --------------------------------------------------
# Load gene sets (Hallmark as default)
# --------------------------------------------------
msig <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

pathways <- split(msig$gene_symbol, msig$gs_name)

# --------------------------------------------------
# Run fgsea
# --------------------------------------------------
gsea_res <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

# --------------------------------------------------
# Tidy and save results
# --------------------------------------------------
gsea_res <- gsea_res %>%
  arrange(padj)

write_tsv(gsea_res, "geomx_GSEA_results.tsv")
