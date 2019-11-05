library(tidyverse)
library(parallel)
source("src/go_enrich.R")

# Permutation test for the enrichment of GO terms associated with the genes
# nearest to significant SNPs versus the background of all genes within the 
# LD-tagged regions. Draw 10k sets of random genes from the 19 regions and 
# test for GO term enrichment.

gene_table <- read_rds("data/gemma/candidate_genes.rds")
nearest <- read_rds("data/gemma/nearest_genes.rds")
go <- read_csv("data/gemma/go_terms.csv")


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
b <- filter(gamer, v2_gene_model %in% bgenes)

cl <- makeCluster(30)
clusterEvalQ(cl, { library(tidyverse); source("src/go_enrich.R") })
clusterExport(cl, list("b"))
enrich <- parApply(cl, gsets, 2, go_enrich)
stopCluster(cl)
write_rds(enrich, "data/gemma/enrich_permute.rds")


enrich <- map(enrich, function(d) filter(d, p_val <= 0.05) %>% pull(Name)) %>%
  unlist(use.names = FALSE) %>%
  tibble(Name = .) %>%
  count(Name) %>%
  mutate(p = n/1e4)

go <- left_join(go, enrich, by = "Name")
