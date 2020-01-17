library(tidyverse)
library(purrrlyr)
library(GenomicRanges)

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
decay_ranges <- with(GM, GRanges(seqnames = Chromosome, 
                                 ranges = IRanges(start = Lower, end = Upper), 
                                 SNP = SNP))
genome_ranges <- with(GM, GRanges(seqnames = Chromosome, 
                                  ranges = IRanges(start = Position, end = Position), 
                                  SNP = SNP))

mirna <- read_delim("~/anno/miRNA.gff3", col_names = FALSE, comment = "#", 
                    delim = "\t") %>%
  select(X1, X4, X5, X9) %>%
  dplyr::rename(chr = X1, start = X4, end = X5, id = X9) %>%
  mutate(chr = str_remove(chr, "Chr") %>% as.integer, 
         id = str_split(id, ";") %>% sapply(., function(x) x[1]) %>%
           str_remove(., "Name="))
mi_ranges <- with(mirna, GRanges(seqnames = chr, 
                                 ranges = IRanges(start = start, end = end), 
                                 ID = id))

lncrna <- read_delim("~/anno/lncRNA.gtf", col_names = FALSE, comment = "#", 
                     delim = "\t") %>%
  select(X1, X4, X5, X9) %>%
  dplyr::rename(chr = X1, start = X4, end = X5, id = X9) %>%
  mutate(chr = str_remove(chr, "chr") %>% as.integer, 
         id = str_split(id, "; ") %>% sapply(., function(x) x[2]) %>%
           str_remove(., "gene_id ") %>% str_remove_all(., '\"')) %>%
  filter(!is.na(chr))
lnc_ranges <- with(lncrna, GRanges(seqnames = chr, 
                                   ranges = IRanges(start = start, end = end), 
                                   ID = id))

cand <- read_rds("data/gemma/nearest_genes.rds")
gene_table <- read_rds("data/gemma/candidate_genes.rds")

snp_ranges <- with(cand, GRanges(seqname = Chromosome, 
                                 ranges = IRanges(start = Position, end = Position), 
                                 ID = SNP))
ld_ranges <- with(gene_table, GRanges(seqname = Chromosome, 
                                      ranges = IRanges(start = Start, end = End), 
                                      Phenotypes = Phenotypes))

# Permutation distributions for overlaps
set.seed(739108)
snps <- sapply(1:1e4, function(i) sample(nrow(GM), nrow(cand), replace = FALSE))
overlaps <- apply(snps, 2, function(s) {
  reduced <- reduce(decay_ranges[s], ignore.strand = TRUE, with.revmap = TRUE)
  mi <- findOverlaps(mi_ranges, reduced, ignore.strand = TRUE)
  lnc <- findOverlaps(lnc_ranges, reduced, ignore.strand = TRUE)
  c(length(unique(mi@from)), 
    length(unique(lnc@from)))
})
over_thresh <- apply(overlaps, 1, quantile, probs = c(0.025, 0.975))

# LD-tagged regions
ld_mi_hit <- findOverlaps(mi_ranges, ld_ranges, ignore.strand = TRUE)
n_mi <- ld_mi_hit@from %>% unique() %>% length() # 14

ld_lnc_hit <- findOverlaps(lnc_ranges, ld_ranges, ignore.strand = TRUE)
n_lnc <- ld_lnc_hit@from %>% unique() %>% length() # 127

# Output tables
mi_table <- mirna[ld_mi_hit@from, ] %>%
  mutate(Phenotypes = ld_ranges@elementMetadata@listData$Phenotypes[ld_mi_hit@to] %>%
           str_replace_all(., ",", ";")) %>%
  dplyr::rename(Chromosome = chr, Start = start, End = end, ID = id) %>%
  arrange(Chromosome, Start, End)
write_csv(mi_table, "data/gemma/STable_miRNA.csv")
lnc_table <- lncrna[ld_lnc_hit@from, ] %>%
  mutate(Phenotypes = ld_ranges@elementMetadata@listData$Phenotypes[ld_lnc_hit@to] %>%
           str_replace_all(., ",", ";")) %>%
  dplyr::rename(Chromosome = chr, Start = start, End = end, ID = id) %>%
  arrange(Chromosome, Start, End)
write_csv(lnc_table, "data/gemma/STable_lncRNA.csv")

# Permutation distributions for distance to nearest
distances <- apply(snps, 2, function(s) {
  cc <- genome_ranges[s]
  mi <- distanceToNearest(cc, mi_ranges, ignore.strand = TRUE)
  lnc <- distanceToNearest(cc, lnc_ranges, ignore.strand = TRUE)
  c(min(mi@elementMetadata@listData$distance), 
    min(lnc@elementMetadata@listData$distance))
})
dist_thresh <- apply(distances, 1, quantile, probs = c(0.025, 0.975))

# SNPs
snp_mi_hit <- findOverlaps(mi_ranges, snp_ranges, ignore.strand = TRUE)   # 0
snp_lnc_hit <- findOverlaps(lnc_ranges, snp_ranges, ignore.strand = TRUE) # 0

snp_mi_near <- distanceToNearest(snp_ranges, mi_ranges, ignore.strand = TRUE)
snp_lnc_near <- distanceToNearest(snp_ranges, lnc_ranges, ignore.strand = TRUE)

sum(snp_mi_near@elementMetadata@listData$distance <= dist_thresh[1, 1])  # 0
sum(snp_mi_near@elementMetadata@listData$distance >= dist_thresh[2, 1])  # 142
sum(snp_lnc_near@elementMetadata@listData$distance <= dist_thresh[1, 2]) # 0
sum(snp_lnc_near@elementMetadata@listData$distance >= dist_thresh[2, 2]) # 142


library(grid)
library(gridExtra)

pA <- ggplot(tibble(Overlaps = overlaps[1, ]), aes(x = Overlaps)) + theme_classic() +
  geom_density() + labs(y = "Density", x = "miRNA LD decay overlaps") + 
  geom_vline(xintercept = over_thresh[1, 1], linetype = 2) +
  geom_vline(xintercept = over_thresh[2, 1], linetype = 2) +
  geom_vline(xintercept = n_mi, linetype = 2, colour = "red")

pB <- ggplot(tibble(Overlaps = overlaps[2, ]), aes(x = Overlaps)) + theme_classic() +
  geom_density() + labs(y = "Density", x = "lncRNA LD decay overlaps") +
  geom_vline(xintercept = over_thresh[1, 2], linetype = 2) +
  geom_vline(xintercept = over_thresh[2, 2], linetype = 2) +
  geom_vline(xintercept = n_lnc, linetype = 2, colour = "red")

pC <- ggplot(tibble(Distance = distances[1, ] + 1), aes(x = Distance)) + theme_classic() + 
  geom_density() + scale_x_log10() +
  labs(y = "Density", x = expression(paste(log[10], "(bp to nearest miRNA)"))) +
  geom_rug(data = tibble(Distance = snp_mi_near@elementMetadata@listData$distance)) +
  geom_vline(aes(xintercept = V1), data = as_tibble(dist_thresh), linetype = 2)

pD <- ggplot(tibble(Distance = distances[2, ] + 1), aes(x = Distance)) + theme_classic() +
  geom_density() + scale_x_log10() + 
  labs(y = "Density", x = expression(paste(log[10], "(bp to nearest lncRNA)"))) +
  geom_rug(data = tibble(Distance = snp_lnc_near@elementMetadata@listData$distance)) + 
  geom_vline(aes(xintercept = V2), data = as_tibble(dist_thresh + 1), linetype = 2)

gA <- grobTree(ggplotGrob(pA), 
               textGrob("A", x = unit(0.03, "npc"), y = unit(0.97, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gB <- grobTree(ggplotGrob(pB), 
               textGrob("B", x = unit(0.03, "npc"), y = unit(0.97, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gC <- grobTree(ggplotGrob(pC), 
               textGrob("C", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gD <- grobTree(ggplotGrob(pD), 
               textGrob("D", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))

lay <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
gp <- arrangeGrob(gA, gB, gC, gD, layout_matrix = lay)
ggsave("figures/select/mi_lnc_enrichment.pdf", plot = gp, width = 10, 
       height = 6, units = "in", dpi = 300)
