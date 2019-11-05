library(tidyverse)
library(purrrlyr)
source("src/go_enrich.R")


# Load data files ---------------------------------------------------------
cand <- read_rds("data/gemma/nearest_genes.rds")
gene_table <- read_rds("data/gemma/candidate_genes.rds") %>%
  by_row(function(r) {
    sapply(unique(cand$Gene), function(g) str_detect(r$Genes, g)) %>%
      sum()
  }, .to = "NGenes") %>%
  unnest(NGenes) %>%
  mutate(Genes = str_split(Genes, ", "))
bgenes <- unlist(gene_table$Genes, use.names = FALSE) %>% unique()
b <- filter(gamer, v2_gene_model %in% bgenes)

# 14 unique combinations of phenotypes -- most have too few entries to be used
# for enrichment analysis. The exception is "lnMSE, TMIN_X0.05_X1.5, TMIN_X0.15_X0.65"
# which has 16 entries.

all_genes <- go_enrich(unique(cand$Gene), b)
all2 <- go_enrich(unique(cand$Gene), gamer)
idx <- sapply(cand$Phenotype, function(p) {
  "lnMSE" %in% p
})
ln_genes <- go_enrich(unique(cand$Gene[idx]), b)

regions <- go_enrich(bgenes, gamer)


all_sig <- filter(all_genes, p_val <= 0.05) %>% pull(Name)
ln_sig <- filter(ln_genes, p_val <= 0.05) %>% pull(Name)
r_sig <- filter(regions, p_val <= 0.05) %>% pull(Name)


gplots::venn(list(Regions = r_sig, Closest = all_sig, MSE = ln_sig))

write_csv(all_genes, "data/gemma/go_terms.csv")
