# GSEA for spatial transcriptomics

library(readr)
library(dplyr)
library(fgsea)
library(msigdbr)

# -----------------------------
# Load DE results
# -----------------------------
de <- read_tsv("spatial_cluster_markers.tsv")

# Rank genes
ranked_genes <- de$avg_log2FC
names(ranked_genes) <- de$gene

ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# -----------------------------
# Load gene sets
# -----------------------------
msig <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

pathways <- split(msig$gene_symbol, msig$gs_name)

# -----------------------------
# Run fgsea
# -----------------------------
fgsea_res <- fgsea(
  pathways = pathways,
  stats = ranked_genes,
  nperm = 10000
)

write_tsv(fgsea_res, "spatial_gsea_results.tsv")
