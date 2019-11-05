library(tidyverse)
library(purrrlyr)
source("src/go_enrich.R")


# Load data files ---------------------------------------------------------
cand <- read_rds("data/single_site_gwas/nearest_genes.rds")$Gene %>%
  unique()


# GO term enrichment ------------------------------------------------------
enrichment <- go_enrich(cand)

write_csv(enrichment, "data/single_site_gwas/go_terms.csv")
