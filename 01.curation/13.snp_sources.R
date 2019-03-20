library(tidyverse)


# Retrieve SNP info for different files -----------------------------------
g2f_snps <- read_rds("data/gbs/recalled_G2F_gbs_snp_info.rds") %>% pull(ID)
ames_snps <- read_rds("data/gbs/recalled_Ames_gbs_snp_info.rds") %>% pull(ID)

rna_snps <- read_tsv("data/gbs/NAM_rna.vcf", comment = "##") %>%
  mutate(`#CHROM` = str_replace(`#CHROM`, "chr", ""), 
         ID = str_replace(ID, "chr", "") %>%
           str_replace(., "-", "_") %>%
           paste0("X", ID)) %>%
  pull(ID)

nam_snps <- read_tsv("data/gbs/NAM/NAM_markers.hmp.txt", comment = "") %>%
  mutate(`rs#` = paste0("X", `rs#`)) %>%
  pull(`rs#`)

hmp1_snps <- read_tsv("data/gbs/HapMap1/maizeHapMapV1.hmp.txt", comment = "") %>%
  mutate(`rs#` = paste0("X", `rs#`)) %>%
  pull(`rs#`)
hmp1_snps <- hmp1_snps[!duplicated(hmp1_snps)]

hmp2_snps <- read_tsv("data/gbs/HapMap2/maizeHapMapV2_recode.vcf", comment = "##") %>%
  mutate(ID = paste0("X", ID)) %>%
  pull(ID)
hmp2_snps <- hmp2_snps[!duplicated(hmp2_snps)]


# SNPs that were retained in the final hybrid table -----------------------
kept <- read_rds("data/gbs/add_snps.rds")$GM %>% pull(SNP)


# Proportion of each source kept ------------------------------------------
tibble(Source = c("G2F", "Ames", "NAM", "HapMap1", "HapMap2", "RNA-seq"), 
       Retained = c(sum(g2f_snps %in% kept)/length(g2f_snps), 
                    sum(ames_snps %in% kept)/length(ames_snps), 
                    sum(nam_snps %in% kept)/length(nam_snps), 
                    sum(hmp1_snps %in% kept)/length(hmp1_snps), 
                    sum(hmp2_snps %in% kept)/length(hmp2_snps), 
                    sum(rna_snps %in% kept)/length(rna_snps)), 
       Removed = 1 - Retained) %>%
  gather(Action, Proportion, Retained:Removed) %>%
  ggplot(., aes(x = Source, y = Proportion, fill = Action)) + theme_classic() +
    geom_col(colour = "black") + labs(x = "", y = "", fill = "") +
    scale_fill_manual(values = c("Retained" = "green", "Removed" = "red")) +
    scale_y_continuous(labels = scales::percent)
ggsave("figures/munge/snp_sources.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Number of sources for each retained SNP ---------------------------------
raw_kept <- c(g2f_snps[g2f_snps %in% kept], ames_snps[ames_snps %in% kept], 
              hmp1_snps[hmp1_snps %in% kept], hmp2_snps[hmp2_snps %in% kept], 
              nam_snps[nam_snps %in% kept], rna_snps[rna_snps %in% kept])

tibble(SNP = raw_kept) %>%
  count(SNP) %>%
  count(n) %>%
  mutate(nn = nn/1e5) %>%
  ggplot(., aes(x = n, y = nn)) + theme_classic() +
    geom_col(colour = "black", fill = "turquoise", alpha = 0.8) +
    labs(x = "Number of Sources", 
         y = expression(paste("Number of SNPs (x", 10^5, ")", sep = ""))) +
    scale_x_continuous(breaks = 1:6, limits = c(0.5, 6.5))
ggsave("figures/munge/source_counts.pdf", width = 6, height = 4, units = "in", dpi = 300)


# For unique SNPs, which source? ------------------------------------------
sources <- tibble(SNP = raw_kept) %>%
  count(SNP) %>%
  filter(n == 1) %>%
  select(-n) %>%
  mutate(Source = if_else(SNP %in% g2f_snps, "G2F", 
                          if_else(SNP %in% ames_snps, "Ames", 
                                  if_else(SNP %in% nam_snps, "NAM", 
                                          if_else(SNP %in% rna_snps, "RNA-seq", 
                                                  if_else(SNP %in% hmp1_snps, "HapMap1", "HapMap2"))))))
count(sources, Source) %>%
  mutate(n = n/sum(n)) %>%
  ggplot(., aes(x = factor(1), y = n)) + theme_classic() +
    geom_col(aes(fill = Source), colour = "black", width = 0.5) +
    labs(x = "", y = "") + scale_fill_brewer(type = "qual", palette = "Dark2") +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() + theme(axis.text.y = element_blank(), 
                         axis.line.y = element_blank(), 
                         axis.ticks.y = element_blank())
ggsave("figures/munge/source_proportion.pdf", width = 6, height = 4, units = "in", dpi = 300)
