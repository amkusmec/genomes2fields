library(tidyverse)
library(mashr)
library(grid)
library(gridExtra)
source("src/manhattan_gemma.R")


m2 <- read_rds("data/gemma/norm_snp_mash.rds")
m3 <- read_rds("data/gemma/norm_snp_mash_all.rds")
phenotypes <- list.files("data/gemma/output", "*\\.assoc\\.txt") %>%
  str_remove(., "norm_") %>%
  str_remove(., "\\.assoc\\.txt")
snps <- read_rds("data/gbs/add_snps.rds")$GM


# Manhattan plot ----------------------------------------------------------
snps <- snps %>% 
  mutate(ps = 1:n(), 
         chr_color = factor(Chromosome %% 2))

bounds <- snps %>%
  group_by(Chromosome) %>%
  summarise(Bound = max(ps)) %>%
  ungroup() %>%
  bind_rows(tibble(Chromosome = 0, Bound = 0), .)

ticks <- snps %>%
  group_by(Chromosome) %>%
  summarise(Tick = (max(ps) - min(ps))/2 + min(ps)) %>%
  ungroup()

assoc <- snps %>%
  filter(SNP %in% rownames(m2$result$lfsr)) %>%
  mutate(sval = apply(m2$result$lfsr, 1, min)) %>%
  filter(sval <= 0.1) %>%
  mutate(sval = -log10(sval))

pA <- ggplot(assoc, aes(x = ps, y = sval, colour = chr_color)) + theme_bw() + 
  geom_vline(aes(xintercept = Bound), bounds, colour = "grey80") + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = -log10(0.1), linetype = 2, colour = "red") +
  scale_colour_manual(values = c("grey60", "black")) + 
  scale_x_continuous(breaks = ticks$Tick, labels = ticks$Chromosome) + 
  labs(colour = "", x = "Chromosome", y = expression(-log[10](s-value))) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "none")


# Phenotypes per SNP ------------------------------------------------------
pB <- get_n_significant_conditions(m2, thresh = 0.1) %>%
  enframe() %>%
  filter(value > 0) %>%
  ggplot(., aes(x = value)) + theme_classic() + 
    geom_histogram(binwidth = 1, colour = "black", fill = "aquamarine", 
                   alpha = 0.8) + 
    labs(x = "# significant phenotypes", y = "# SNPs")


# Dominant covariance component -------------------------------------------
# ED1 (relative frequency = 36.2%); ED2 (relative frequency = 11.7%)
pi_mat <- matrix(m2$fitted_g$pi[-1], ncol = length(m2$fitted_g$Ulist), 
                 nrow = length(m2$fitted_g$grid), byrow = TRUE)
colnames(pi_mat) <- names(m2$fitted_g$Ulist)
colSums(pi_mat) %>% sort(decreasing = TRUE) %>% head()

cmat <- cov2cor(m2$fitted_g$Ulist$ED_1)
rownames(cmat) <- colnames(cmat) <- colnames(m2$result$lfsr)
cmat[upper.tri(cmat, diag = FALSE)] <- NA

pC <- as_tibble(cmat, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = unique(Phenotype1), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(unique(Phenotype2)), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() + 
    geom_tile() +  scale_x_discrete(position = "top") + 
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625)) + 
    labs(x = "", y = "", fill = "", subtitle = "relative frequency = 36.2%")


# First eigenvector -------------------------------------------------------
vold <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))$v
u <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))$u
d <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))$d
v <- vold[order(colnames(m2$result$lfsr)), ]

pD <- tibble(Phenotype = colnames(m2$result$lfsr), 
             V = v[, 1]/v[, 1][which.max(abs(v[, 1]))]) %>%
  mutate(Variable = str_remove(Phenotype, "_X[01]\\.[0-9]{2,3}_X[01]\\.[0-9]{1,3}")) %>%
  ggplot(., aes(x = Phenotype, y = V)) + theme_classic() + 
    geom_col(aes(fill = Variable), colour = "black") + 
    labs(x = "", y = "", fill = "", 
         subtitle = paste0("Eigenvector 1, (PVE = ", 
                           round(d[1]^2/sum(d^2), 2)*100, "%)")) + 
    scale_fill_manual(values = c("Hybrid" = "grey60", "TMAX" = "red", 
                                 "SR" = "orange", "NET" = "brown")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Sharing by sign ---------------------------------------------------------
m_data <- read_rds("data/gemma/mash_data.rds")
m_data[-3] <- lapply(m_data[-3], function(m) m[rownames(m2$result$lfsr), ])
pm_mash_beta <- m2$result$PosteriorMean*m_data$s_hat
pm_mash_beta_mag <- pm_mash_beta[rowSums(m2$result$lfsr < 0.1) > 0, ]
lfsr_mash <- m2$result$lfsr[rowSums(m2$result$lfsr < 0.1) > 0, ]
shared_fold_size <- matrix(NA, nrow = ncol(lfsr_mash), ncol = ncol(lfsr_mash))
colnames(shared_fold_size) <- rownames(shared_fold_size) <- colnames(m2$result$lfsr)

for (i in 1:ncol(lfsr_mash)) {
  for (j in 1:ncol(lfsr_mash)) { 
    sig_row <- which(lfsr_mash[, i] < 0.1)
    sig_col <- which(lfsr_mash[, j] < 0.1)
    a <- union(sig_row, sig_col)
    quotient <- pm_mash_beta_mag[a, i]/pm_mash_beta_mag[a, j]
    shared_fold_size[i, j] <- mean(quotient > 0)
  }
}

shared_fold_size <- shared_fold_size[order(rownames(shared_fold_size)), 
                                     order(colnames(shared_fold_size))]
shared_fold_size[upper.tri(shared_fold_size, diag = FALSE)] <- NA

pE <- as_tibble(shared_fold_size, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() + 
    geom_tile() + labs(x = "", y = "", fill = "") + 
    scale_fill_distiller(type = "seq", palette = "Greens", 
                         limits = c(0, 1), direction = 1, labels = scales::percent) + 
    scale_x_discrete(position = "top") + 
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))


# Assemble the figure -----------------------------------------------------
gA <- grobTree(ggplotGrob(pA), 
               textGrob("A", x = unit(0.01, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gB <- grobTree(ggplotGrob(pB), 
               textGrob("E", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gC <- grobTree(ggplotGrob(pC), 
               textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gD <- grobTree(ggplotGrob(pD), 
               textGrob("C", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gE <- grobTree(ggplotGrob(pE), 
               textGrob("D", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))

lay <- matrix(c(1, 1, 2, 4, 3, 5), ncol = 2, byrow = TRUE)
gp <- arrangeGrob(gA, gC, gD, gE, gB, layout_matrix = lay)
ggsave("figures/select/mash_plots.pdf", plot = gp, width = 12, height = 12, 
       units = "in", dpi = 300)
