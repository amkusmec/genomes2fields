library(tidyverse)
library(GenomicRanges)
library(parallel)
source("src/go_enrich.R")

decay <- read_csv("data/decay_regions.csv")
GM <- read_csv("data/gbs/add_map.csv") %>%
  by_row(function(r) {
    i <- which(r$Chromosome[1] == decay$Chr &
                 r$Position[1] >= decay$Start &
                 r$Position[1] < decay$End)
    decay$Decay[i]
  }, .to = "Decay") %>%
  unnest() %>%
  mutate(Lower = Position - Decay,
         Lower = if_else(Lower < 0, 0, Lower),
         Upper = Position + Decay)
snpRanges <- with(GM, GRanges(seqnames = Chromosome, 
                              ranges = IRanges(start = Lower, end = Upper), 
                              SNP = SNP))

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

set.seed(375098)
snps <- sapply(1:1e4, function(i) sample(nrow(GM), 44, replace = FALSE))
genes <- apply(snps, 2, function(s) {
  reduced <- reduce(snpRanges[s], ignore.strand = TRUE, with.revmap = TRUE)
  genes <- findOverlaps(geneRanges, reduced, ignore.strand = TRUE)
  geneRanges$ids[unique(genes@from)]
})

cl <- makeCluster(10)
clusterEvalQ(cl, source("src/go_enrich.R"))

e <- parLapply(cl, genes, function(g) {
  enrich <- go_enrich(g)
  filter(enrich, q_val <= 0.05) %>% pull(Name)
})

stopCluster(cl)

write_rds(e, "data/go_permute.rds")
