library(tidyverse)
library(mgcv)
library(QGenTools)
library(purrrlyr)
library(GenomicRanges)

est_distance <- function(ld, threshold = 0.3, k = 30) {
  require(mgcv)
  
  init_gam <- gam(R2 ~ s(Bin, k = k, bs = "cr"), data = ld)
  sm <- smoothCon(s(Bin, k = k, bs = "cr"), ld, knots = NULL)[[1]]
  mc <- mono.con(sm$xp, up = FALSE)
  M <- list(X = sm$X, y = ld$R2, C = matrix(0, 0, 0), Ain = mc$A, 
            bin = mc$b, sp = init_gam$sp, p = -sm$xp, S = sm$S, 
            w = ld$R2*0 + 1, off = 0)
  p <- pcls(M)
  if (all(is.nan(p))) {
    min(ld$Bin)
  } else {
    x <- seq(min(ld$Bin), max(ld$Bin), 1e3)
    y <- abs(drop(Predict.matrix(sm, data.frame(Bin = x)) %*% p) - threshold)
    x[which.min(y)]
  }
}


m2 <- read_rds("data/gemma/norm_snp_mash.rds")
sig <- which(rowSums(apply(m2$result$lfsr, 2, function(x) x <= 0.05)) > 0)
sig <- tibble(SNP = rownames(m2$result$lfsr[sig, ])) %>%
  mutate(Chromosome = str_remove(SNP, "X") %>%
           str_remove(., "_[0-9]*") %>% as.integer(), 
         Position = str_remove(SNP, "X[0-9]{1,2}_") %>%
           as.integer(), 
         Up = NA, Down = NA) %>%
  arrange(Chromosome, Position)

snps <- read_rds("data/gbs/add_snps.rds")
AB <- realized_ab(snps$GD)
struct <- read_rds("data/gbs/pca_covariates.rds")
window <- 500

for (s in 1:nrow(sig)) {
  cat(s, "/", nrow(sig), "\r")
  
  # Find the location of the SNP in the genome
  chr <- which(snps$GM$Chromosome == sig$Chromosome[s])
  idx <- which(snps$GM$SNP == sig$SNP[s])
  idx1 <- if_else(idx - window < min(chr), min(chr), as.integer(idx - window))
  idx2 <- if_else(idx + window < max(chr), as.integer(idx + window), max(chr))
  
  # LD decay upstream of the SNP
  r2_up <- r2vs(snps$GD[, idx1:idx], AB, struct) %>%
    mutate_at(c("Locus1", "Locus2"), as.character) %>%
    filter(Locus1 == sig$SNP[s] | Locus2 == sig$SNP[s]) %>%
    mutate(Locus1 = str_remove(Locus1, "X[0-9]{1,2}_") %>% as.integer, 
           Locus2 = str_remove(Locus2, "X[0-9]{1,2}_") %>% as.integer, 
           Bin = abs(Locus1 - Locus2)) %>%
    arrange(Bin)
  sig$Up[s] <- est_distance(r2_up, threshold = 0.05)
  
  # LD decay downstream of the SNP
  r2_down <- r2vs(snps$GD[, idx:idx2], AB, struct) %>%
    mutate_at(c("Locus1", "Locus2"), as.character) %>%
    filter(Locus1 == sig$SNP[s] | Locus2 == sig$SNP[s]) %>%
    mutate(Locus1 = str_remove(Locus1, "X[0-9]{1,2}_") %>% as.integer,
           Locus2 = str_remove(Locus2, "X[0-9]{1,2}_") %>% as.integer, 
           Bin = abs(Locus1 - Locus2)) %>%
    arrange(Bin)
  sig$Down[s] <- est_distance(r2_down, threshold = 0.05)
}

sig <- sig %>%
  by_row(function(r) {
    s <- which(m2$result$lfsr[r$SNP[1], ] <= 0.05)
    colnames(m2$result$lfsr)[s]
  }, .to = "Phenotype")
sig$N_pheno <- sapply(sig$Phenotype, length)


sranges <- with(sig, GRanges(seqnames = Chromosome, 
                             ranges = IRanges(start = Position - Up, end = Position + Down), 
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

# 36/38 (94.7%) ranges contain genes
genes <- findOverlaps(geneRanges, reduced, ignore.strand = TRUE)
idx1 <- genes@from; idx2 <- genes@to

