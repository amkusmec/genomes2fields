library(tidyverse)
library(parallel)
source("src/go_enrich_at.R")


ld <- read_csv("data/decay_regions.csv")
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
gene_ranges <- with(fgs, GRanges(seqnames = chr, 
                                 ranges = IRanges(start = start, end = end), 
                                 GID = gene))

set.seed(849671)
gsets <- lapply(1:1e4, function(i) {
  idx <- sample(nrow(snps), 144)
  temp <- snps[idx, ] %>%
    by_row(function(r) {
      i <- which(r$Chromosome[1] == ld$Chr &
                   r$Position[1] >= ld$Start & 
                   r$Position[1] < ld$End)
      ld$Decay[i]
    }, .to = "Decay") %>%
    unnest() %>%
    mutate(Lower = Position - Decay, 
           Upper = Position + Decay)
  snp_ranges <- with(temp, GRanges(seqnames = Chromosome, 
                                   ranges = IRanges(start = Lower, end = Upper)))
  reduced <- reduce(snp_ranges)
  hits <- findOverlaps(gene_ranges, reduced, ignore.strand = TRUE)
  gene_ranges@elementMetadata@listData$GID[hits@from] %>% unique()
})


cl <- makeCluster(30)
clusterEvalQ(cl, { library(tidyverse); library(purrrlyr); source("src/go_enrich_at.R") })
enrich <- parLapply(cl, gsets, go_enrich)
stopCluster(cl)
write_rds(enrich, "data/gemma/ld_vs_genome_permute.rds")
