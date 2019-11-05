library(tidyverse)
source("src/go_enrich.R")

term_counts <- gamer %>%
  count(GO) %>%
  inner_join(., obo, by = "GO")


perm <- read_rds("data/go_permute.rds")
sum(sapply(perm, length) == 0)

non_zero <- perm[sapply(perm, length) > 0]
summary(sapply(non_zero, length))
quantile(sapply(non_zero, length), probs = seq(0.1, 0.9, 0.1))
terms <- unique(unlist(non_zero, use.names = FALSE))
rand_term_counts <- filter(term_counts, Name %in% terms)
rand_term_counts <- unlist(non_zero, use.names = FALSE) %>%
  table() %>% enframe() %>%
  inner_join(., rand_term_counts, by = c("name" = "Name")) %>%
  mutate(value = as.integer(value))

library(GenomicRanges)
gff <- read_delim("~/anno/ZmB73_5b_FGS.gff", comment = "#", delim = "\t", 
                  progress = FALSE, col_names = FALSE) %>%
  filter(X3 == "gene", !is.na(X1)) %>%
  select(X1, X4, X5, X9) %>%
  mutate(X9 = str_split(X9, ";") %>% sapply(., function(x) x[1]), 
         X9 = str_remove(X9, "ID=")) %>%
  arrange(X1, X4) %>%
  dplyr::rename(chr = X1, start = X4, end = X5, gene = X9)
geneRanges <- with(gff, GRanges(seqname = chr, 
                                ranges = IRanges(start = start, end = end), 
                                ids = gene))

hormone_genes <- unique(gamer$v2_gene_model[gamer$GO == "GO:0060416"])
hormone_ranges <- geneRanges[geneRanges$ids %in% hormone_genes]

water_genes <- unique(gamer$v2_gene_model[gamer$GO == "GO:0080148"])
water_ranges <- geneRanges[geneRanges$ids %in% water_genes]
