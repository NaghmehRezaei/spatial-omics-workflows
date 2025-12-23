# Generate PI-facing summary tables for GeoMx spatial transcriptomics analysis
# Inputs: limma DE results and GSEA results
# Outputs: concise, ranked tables for interpretation and reporting

library(readr)
library(dplyr)

# --------------------------------------------------
# Load differential expression results
# --------------------------------------------------
de <- read_tsv("geomx_DE_limma_results.tsv")

# --------------------------------------------------
# Create DE summary table
# --------------------------------------------------
de_summary <- de %>%
  mutate(
    Direction = ifelse(logFC > 0, "Upregulated", "Downregulated")
  ) %>%
  arrange(adj.P.Val) %>%
  select(
    Gene,
    logFC,
    AveExpr,
    t,
    P.Value,
    Adjusted_P_value = adj.P.Val,
    Direction
  ) %>%
  slice(1:50)

write_tsv(de_summary, "PI_DE_summary_top50.tsv")

# --------------------------------------------------
# Load GSEA results
# --------------------------------------------------
gsea <- read_tsv("geomx_GSEA_results.tsv")

# --------------------------------------------------
# Create pathway summary table
# --------------------------------------------------
pathway_summary <- gsea %>%
  mutate(
    Direction = ifelse(NES > 0, "Enriched_in_Group1", "Enriched_in_Group2")
  ) %>%
  arrange(padj) %>%
  select(
    Pathway = pathway,
    NES,
    Direction,
    Adjusted_P_value = padj,
    Pathway_Size = size
  ) %>%
  slice(1:25)

write_tsv(pathway_summary, "PI_pathway_summary_top25.tsv")
