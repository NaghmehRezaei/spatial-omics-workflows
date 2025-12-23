# Generate PI-facing summary tables for spatial analysis

library(readr)
library(dplyr)

gsea <- read_tsv("spatial_gsea_results.tsv")

summary_tbl <- gsea %>%
  mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
  arrange(padj) %>%
  select(
    Pathway = pathway,
    NES,
    Direction,
    Adjusted_P_value = padj,
    Pathway_Size = size
  ) %>%
  slice(1:25)

write_tsv(summary_tbl, "PI_spatial_pathway_summary.tsv")
