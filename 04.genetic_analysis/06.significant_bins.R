### Identify candidate genes based on chromosome-wise LD decay.

library(tidyverse)
library(readxl)
library(QGenTools)
library(mgcv)
library(purrrlyr)
library(GenomicRanges)
library(parallel)


map <- read_xlsx("data/nam_recomb_bins.xlsx", skip = 1) %>%
  group_by(Chr) %>%
  mutate(Gen = c(`Genetic map (cM)`[1], 
                 diff(`Genetic map (cM)`, differences = 1))) %>%
  ungroup() %>%
  mutate(R = Gen/`Bin size (Mb)`) %>%
  split(., .$Chr)

breaks <- map_df(map, function(x) {
  cutoff <- 1:floor(quantile(x$R, probs = 0.96))
  runs <- list()
  for (xx in cutoff) {
    idx <- which(x$R <= xx)
    r <- rle(diff(idx))
    rmax <- which.max(r$lengths)
    rlen <- sum(r$lengths[1:(rmax - 1)])
    runs[[xx]] <- idx[(rlen + 1):(rlen + 1 + r$lengths[rmax])]
  }
  
  idx <- which.max(sapply(runs, length))
  tibble(Chr = x$Chr[1],
         Start = x$Start[min(runs[[idx]])],
         End = x$Stop[max(runs[[idx]])])
}) %>%
  mutate(Start = Start*1e6, 
         End = End*1e6)

breaks$End[5] <- map[[5]]$Stop[352]*1e6
breaks$End[10] <- map[[10]]$Stop[254]*1e6

centro <- read_tsv("~/anno/Centro.gff3", col_names = FALSE, comment = "##") %>%
  mutate(X1 = str_remove(X1, "Chr") %>% as.integer(), 
         X4 = X4/1e6, 
         X5 = X5/1e6) %>%
  dplyr::rename(Chr = X1)

label_chr <- function(x) paste("Chromosome ", x)

bind_rows(map) %>%
  ggplot(., aes(x = Start, y = R)) + theme_bw() + geom_step() +
    facet_wrap(~ Chr, scales = "free", ncol = 5, 
               labeller = labeller(Chr = label_chr)) + 
    geom_vline(aes(xintercept = Start), mutate(breaks, Start = Start/1e6), 
               linetype = 2, colour = "red") +
    geom_vline(aes(xintercept = End), mutate(breaks, End = End/1e6), 
               linetype = 2, colour = "red") + 
    geom_vline(aes(xintercept = X4), centro, linetype = 2, colour = "green") + 
    geom_vline(aes(xintercept = X5), centro, linetype = 2, colour = "green") + 
    labs(x = "Physical Position (Mb)", y = "cM/Mb")
ggsave("figures/chrom_split.pdf", width = 10, height = 5, units = "in", dpi = 300)

m2 <- read_rds("data/gemma/norm_snp_mash.rds")
sig <- which(rowSums(apply(m2$result$lfsr, 2, function(x) x <= 0.1)) > 0)
sig <- tibble(SNP = rownames(m2$result$lfsr[sig, ])) %>%
  mutate(Chromosome = str_remove(SNP, "X") %>%
           str_remove(., "_[0-9]*") %>% as.integer(), 
         Position = str_remove(SNP, "X[0-9]{1,2}_") %>% as.integer()) %>%
  arrange(Chromosome, Position)

snps <- read_rds("data/gbs/add_snps.rds")
AB <- sommer::A.mat(snps$GD - 1)
struct <- read_rds("data/gbs/pca_covariates.rds")

bin_size <- 1e3
window_size <- 500
step_size <- 50

cl <- makeCluster(10)
clusterEvalQ(cl, { library(tidyverse); library(QGenTools); library(mgcv) })
clusterExport(cl, list("snps", "AB", "struct", "bin_size", "window_size", 
                       "step_size", "breaks"))

ld <- parLapply(cl, 1:10, function(r) {
  ld_temp <- tibble(Chr = r, 
                    Start = c(1, breaks$Start[r], breaks$End[r]), 
                    End = c(breaks$Start[r], breaks$End[r], 
                            max(snps$GM$Position[snps$GM$Chromosome == r]) + 1), 
                    Decay = 0)
  
  for (i in 1:nrow(ld_temp)) {
    chr <- which(snps$GM$Chromosome == r & snps$GM$Position >= ld_temp$Start[i] &
                   snps$GM$Position < ld_temp$End[i])
    r2 <- list()
    counter <- 1
    start <- 1
    end <- start + window_size - 1
    
    repeat {
      r2[[counter]] <- r2vs(snps$GD[, chr[start:end]], AB, struct) %>%
        mutate_at(c("Locus1", "Locus2"), as.character) %>%
        mutate(Locus1 = str_remove(Locus1, "X[0-9]{1,2}_") %>% as.integer,
               Locus2 = str_remove(Locus2, "X[0-9]{1,2}_") %>% as.integer, 
               Bin = floor(abs(Locus1 - Locus2)/bin_size))
      
      if (end + step_size >= length(chr)) break
      counter <- counter + 1
      start <- start + step_size
      end <- end + step_size
    }
    
    r2 <- do.call("rbind", r2) %>%
      group_by(Bin) %>%
      summarise(R2 = mean(R2)) %>%
      ungroup()
    
    init_gam <- gam(R2 ~ s(Bin, k = 100, bs = "cr"), data = r2)
    y <- drop(predict.gam(init_gam))
    ld_temp$Decay[i] <- r2$Bin[which(y <= 0.05)[1]]*bin_size
  }
  
  return(ld_temp)
})

ld <- bind_rows(ld) %>%
  mutate(Region = rep(c("Short arm", "Pericentromeric region", "Long arm"), times = 10))
stopCluster(cl)
write_csv(ld, "data/decay_regions.csv")

ld %>%
  mutate(Chr = factor(Chr), 
         Decay = Decay/1e6,
         Region = factor(Region, levels = c("Short arm", "Pericentromeric region", "Long arm"), ordered = TRUE)) %>%
ggplot(., aes(x = Chr, y = Decay, group = Region, fill = Region)) + 
  geom_col(position = "dodge", colour = "black") + theme_bw() + #scale_y_log10() +
  #labs(x = "Chromosome", y = expression(paste(log[10], "(Decay Distance)"))) +
  labs(x = "Chromosome", y = "Decay Distance (Mb)", fill = "") +
  scale_fill_brewer(type = "qual", palette = "Set2")
ggsave("figures/select/ld_decay.pdf", width = 6, height = 4, units = "in", dpi = 300)

sig <- by_row(sig, function(r) {
  i <- which(r$Chromosome[1] == ld$Chr &
               r$Position[1] >= ld$Start &
               r$Position[1] < ld$End)
  ld$Decay[i]
}, .to = "Decay") %>% 
  select(SNP:N_pheno) %>% 
  # unnest(cols = c(Decay)) %>%
  mutate(Lower = Position - Decay, 
         Upper = Position + Decay) %>%
  # by_row(function(r) {
  #   s <- which(m2$result$lfsr[r$SNP[1], ] <= 0.1)
  #   colnames(m2$result$lfsr)[s]
  # }, .to = "Phenotype2") %>%
  mutate(N_pheno = sapply(Phenotype, length))

sig %>%
  select(SNP:Position, Lower:Upper, N_pheno, Phenotype) %>%
  mutate(Phenotype = sapply(Phenotype, paste, collapse = "; ")) %>%
  dplyr::rename(`# Phenotypes` = N_pheno) %>% 
  write_csv("data/gemma/significant_SNPs.csv")

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
idx1 <- genes@from; length(unique(idx1)) # 3,458 genes
idx2 <- genes@to; length(unique(idx2))   # 28/28 (100%) ranges contain genes

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
write_rds(gene_table, "data/gemma/candidate_genes.rds")

gene_table %>%
  mutate(Genes = str_split(Genes, ", ")) %>%
  unnest(Genes) %>%
  dplyr::rename(Gene = Genes) %>%
  inner_join(., gff, by = c("Chromosome" = "chr", "Gene" = "gene")) %>%
  select(Gene, Chromosome, start, end, SNPs, Phenotypes) %>%
  dplyr::rename(Start = start, End = end) %>%
  mutate(SNPs = str_replace_all(SNPs, ",", ";"), 
         Phenotypes = str_replace_all(Phenotypes, ",", ";")) %>%
  write_csv("data/gemma/candidate_genes.csv")


snp_ranges <- with(sig, GRanges(seqnames = Chromosome, 
                                ranges = IRanges(start = Position, end = Position), 
                                SNP = SNP))
dist_to_nearest <- distanceToNearest(snp_ranges, geneRanges, ignore.strand = TRUE)
sig$Gene[dist_to_nearest@from] <- geneRanges$ids[dist_to_nearest@to]
sig$Distance <- dist_to_nearest@elementMetadata@listData$distance
sum(sig$Distance == 0) # 26/144 (18.1%) of SNPs are in a gene
length(unique(sig$Gene)) # 55 unique genes
max(sig$Distance) # Farthest gene is 309,713 bp from a SNP

write_rds(sig, "data/gemma/nearest_genes.rds")

sig %>%
  by_row(function(r) {
    paste(r$Phenotype[[1]], collapse = "; ")
  }, .to = "P", .collate = "rows") %>%
  dplyr::select(SNP, Chromosome, Position, Gene, Distance, N_pheno, P) %>%
  rename(`Distance to SNP (bp)` = Distance, `# Phenotypes` = N_pheno, 
         Phenotypes = P) %>%
  write_csv(., "data/gemma/nearest_genes.csv")

