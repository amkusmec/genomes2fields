library(tidyverse)
library(purrrlyr)
source("src/go_enrich_at.R")


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
b <- filter(anno, GID %in% bgenes)

sapply(cand$Phenotype, function(l) paste(l, collapse = ",")) %>% unique()
# 9 unique combinations of phenotypes
sapply(cand$Phenotype, function(l) paste(l, collapse = ",")) %>% table()

closest_vs_ld <- go_enrich(unique(cand$Gene), b)
closest_vs_genome <- go_enrich(unique(cand$Gene), anno)
ld_vs_genome <- go_enrich(bgenes, anno)


cl_sig <- filter(closest_vs_ld, p_val <= 0.05) %>% pull(Name)
cg_sig <- filter(closest_vs_genome, p_val <= 0.05) %>% pull(Name)
lg_sig <- filter(ld_vs_genome, p_val <= 0.05) %>% pull(Name)


gplots::venn(list(CL = cl_sig, CG = cg_sig, LG = lg_sig))

write_csv(closest_vs_ld, "data/gemma/closest_vs_ld_go_terms.csv")
write_csv(closest_vs_genome, "data/gemma/closest_vs_genome_go_terms.csv")
write_csv(ld_vs_genome, "data/gemma/ld_vs_genome_go_terms.csv")
