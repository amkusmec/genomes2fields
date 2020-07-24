### Permutation analysis for enrichment of GO terms at genes closest to the 
### significant SNPs against the background of all genes.

library(tidyverse)
library(parallel)
source("src/go_enrich_at.R")


nearest <- read_rds("data/gemma/nearest_genes.rds")
fgs <- read_delim("~/anno/ZmB73_5b_FGS.gff", comment = "#", delim = "\t", 
                  progress = FALSE, col_names = FALSE) %>%
  filter(X3 == "gene", !is.na(X1)) %>%
  select(X1, X4, X5, X9) %>%
  mutate(X9 = str_split(X9, ";") %>% sapply(., function(x) x[1]), 
         X9 = str_remove(X9, "ID=")) %>%
  arrange(X1, X4) %>%
  rename(chr = X1, start = X4, end = X5, gene = X9)
snps <- read_rds("data/gbs/add_snps.rds")$GM

library(GenomicRanges)
snp_ranges <- with(snps, GRanges(seqnames = Chromosome, 
                                 ranges = IRanges(start = Position, end = Position)))
gene_ranges <- with(fgs, GRanges(seqnames = chr, 
                                 ranges = IRanges(start = start, end = end), 
                                 GID = gene))

set.seed(328075)
gsets <- lapply(1:1e4, function(i) {
  idx <- nearest(snp_ranges[sample(length(snp_ranges), 144)], gene_ranges, 
                 select = "all", ignore.strand = TRUE)
  gene_ranges@elementMetadata@listData$GID[idx@to] %>% unique()
})


cl <- makeCluster(30)
clusterEvalQ(cl, { library(tidyverse); library(purrrlyr); source("src/go_enrich_at.R") })
enrich <- parLapply(cl, gsets, go_enrich)
stopCluster(cl)
write_rds(enrich, "data/gemma/closest_vs_genome_permute.rds")
