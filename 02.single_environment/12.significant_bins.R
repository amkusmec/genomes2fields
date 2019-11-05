library(tidyverse)
library(QGenTools)
library(purrrlyr)
library(GenomicRanges)

m2 <- read_rds("data/single_site_gwas/common_snp_mash.rds")
sig <- which(rowSums(apply(m2$result$lfsr, 2, function(x) x <= 0.05)) > 0)
sig <- tibble(SNP = rownames(m2$result$lfsr[sig, ])) %>%
  mutate(Chromosome = str_remove(SNP, "X") %>%
           str_remove(., "_[0-9]*") %>% as.integer(), 
         Position = str_remove(SNP, "X[0-9]{1,2}_") %>%
           as.integer()) %>%
  arrange(Chromosome, Position)

# snps <- read_rds("data/gbs/add_snps.rds")
# AB <- sommer::A.mat(snps$GD - 1)
# struct <- read_rds("data/gbs/pca_covariates.rds")
# window <- 500
# 
# for (s in 1:nrow(sig)) {
#   cat(s, "/", nrow(sig), "\r")
#   
#   # Find the location of the SNP in teh genome
#   chr <- which(snps$GM$Chromosome == sig$Chromosome[s])
#   idx <- which(snps$GM$SNP == sig$SNP[s])
#   idx1 <- if_else(idx - window < min(chr), min(chr), as.integer(idx - window))
#   idx2 <- if_else(idx + window < max(chr), as.integer(idx + window), max(chr))
#   
#   # LD decay upstream of the SNP
#   r2_up <- r2vs(snps$GD[, idx1:idx], AB, struct) %>%
#     mutate_at(c("Locus1", "Locus2"), as.character) %>%
#     filter(Locus1 == sig$SNP[s] | Locus2 == sig$SNP[s]) %>%
#     mutate(Locus1 = str_remove(Locus1, "X[0-9]{1,2}_") %>% as.integer, 
#            Locus2 = str_remove(Locus2, "X[0-9]{1,2}_") %>% as.integer, 
#            Bin = abs(Locus1 - Locus2)) %>%
#     arrange(Bin)
#   sig$Up[s] <- est_distance(r2_up, threshold = 0.05, k = 30)
#   
#   # LD decay downstream of the SNP
#   r2_down <- r2vs(snps$GD[, idx:idx2], AB, struct) %>%
#     mutate_at(c("Locus1", "Locus2"), as.character) %>%
#     filter(Locus1 == sig$SNP[s] | Locus2 == sig$SNP[s]) %>%
#     mutate(Locus1 = str_remove(Locus1, "X[0-9]{1,2}_") %>% as.integer, 
#            Locus2 = str_remove(Locus2, "X[0-9]{1,2}_") %>% as.integer, 
#            Bin = abs(Locus1 - Locus2)) %>%
#     arrange(Bin)
#   sig$Down[s] <- est_distance(r2_down, threshold = 0.05, k = 30)
# }

ld <- read_csv("data/decay_regions.csv")
sig <- by_row(sig, function(r) {
  i <- which(r$Chromosome[1] == ld$Chr &
               r$Position[1] >= ld$Start &
               r$Position[1] < ld$End)
  ld$Decay[i]
}, .to = "Decay") %>%
  unnest() %>%
  mutate(Lower = Position - Decay, 
         Upper = Position + Decay) %>%
  by_row(function(r) {
    s <- which(m2$result$lfsr[r$SNP[1], ] <= 0.05)
    colnames(m2$result$lfsr)[s]
  }, .to = "Phenotype") %>%
  mutate(N_pheno <- sapply(Phenotype, length))

sranges <- with(sig, GRanges(seqnames = Chromosome, 
                             ranges = IRanges(start = Lower, end = Upper), 
                             SNP = SNP))
reduced <- reduce(sranges, ignore.strand = TRUE, with.revmap = TRUE)

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

genes <- findOverlaps(geneRanges, reduced, ignore.strand = TRUE)
idx1 <- genes@from; length(unique(idx1)) # 23,830 genes
idx2 <- genes@to; length(unique(idx2)) # 119/119 (100%) ranges contain genes

gene_table <- as.data.frame(reduced) %>% as_tibble() %>%
  mutate(seqnames = as.character(seqnames) %>% as.integer(), 
         SNPs = NA, Genes = NA, Phenotypes = NA) %>%
  select(-width, -strand)

for (i in 1:nrow(gene_table)) {
  s <- eval(gene_table$revmap[[i]]) %>% sort()
  gene_table$SNPs[i] <- paste(sig$SNP[s], collapse = ", ")
  gene_table$Genes[i] <- geneRanges$ids[idx1[which(idx2 == i)]] %>%
    paste(., collapse = ", ")
  gene_table$Phenotypes[i] <- sig$Phenotype[s] %>% unlist() %>%
    unique() %>% sort() %>% paste(., collapse = ", ")
}

gene_table <- gene_table %>%
  select(-revmap) %>%
  dplyr::rename(Chromosome = seqnames, Start = start, End = end)
write_rds(gene_table, "data/single_site_gwas/candidate_genes.rds")

snp_ranges <- with(sig, GRanges(seqnames = Chromosome, 
                                ranges = IRanges(start = Position, end = Position), 
                                SNP = SNP))
dist_to_nearest <- distanceToNearest(snp_ranges, geneRanges, ignore.strand = TRUE)
sig$Gene[dist_to_nearest@from] <- geneRanges$ids[dist_to_nearest@to]
sig$Distance <- dist_to_nearest@elementMetadata@listData$distance

sum(sig$Distance == 0)/nrow(sig) # 1,220/3,377 (36%) of SNPs are in a gene
length(unique(sig$Gene)) # 1,418 unique genes
max(sig$Distance) # Farthest gene is ~165kb from a SNP

write_rds(sig, "data/single_site_gwas/nearest_genes.rds")
