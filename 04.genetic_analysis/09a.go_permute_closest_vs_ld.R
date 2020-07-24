library(tidyverse)
library(parallel)
source("src/go_enrich_at.R")


gene_table <- read_rds("data/gemma/candidate_genes.rds")
nearest <- read_rds("data/gemma/nearest_genes.rds")
closest_vs_ld <- read_csv("data/gemma/closest_vs_ld_go_terms.csv")

gene_table <- gene_table %>%
  by_row(function(r) {
    sapply(unique(nearest$Gene), function(g) str_detect(r$Genes, g)) %>%
      sum()
  }, .to = "NGenes") %>%
  unnest(NGenes) %>%
  mutate(Genes = str_split(Genes, ", "))

set.seed(772189)
gsets <- sapply(1:1e4, function(i) {
  map2(gene_table$Genes, gene_table$NGenes, function(x, y) {
    sample(x, y, replace = FALSE)
  }) %>% unlist()
})

bgenes <- unlist(gene_table$Genes, use.names = FALSE) %>% unique()
b <- filter(anno, GID %in% bgenes)

cl <- makeCluster(30)
clusterEvalQ(cl, { library(tidyverse); library(purrrlyr); source("src/go_enrich_at.R") })
clusterExport(cl, list("b"))
enrich <- parApply(cl, gsets, 2, go_enrich)
stopCluster(cl)
write_rds(enrich, "data/gemma/closest_vs_ld_permute.rds")
