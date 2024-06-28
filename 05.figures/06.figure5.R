library(tidyverse)
library(mashr)
library(grid)
library(gridExtra)
source("src/manhattan_gemma.R")


m2 <- read_rds("data/gemma/norm_snp_mash.rds")
m3 <- read_rds("data/gemma/norm_snp_mash_all.rds")
phenotypes <-c("Intercept", "Whole season net\nevapotranspiration", 
               "Mid-season\nsolar radiation", 
               "Pre-anthesis\nmax. temp.", "Post-anthesis\nmax. temp.")
fills <- c("Intercept" = "grey60", 
           "Whole season net\nevapotranspiration" = "brown", 
           "Mid-season\nsolar radiation" = "orange", 
           "Pre-anthesis\nmax. temp." = "red", 
           "Post-anthesis\nmax. temp." = "red")
fills2 <- fills[order(names(fills))]
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

sig <- which(rowSums(apply(m2$result$lfsr, 2, function(x) x <= 0.1)) > 0)
assoc <- tibble(SNP = rownames(m2$result$lfsr[sig, ])) %>%
  mutate(Chromosome = str_remove(SNP, "X") %>%
           str_remove("_[0-9]*") %>% as.integer(), 
         Position = str_remove(SNP, "X[0-9]{1,2}_") %>% as.integer(), 
         sval = apply(m2$result$lfsr[sig, ], 1, min)) %>% 
  filter(sval <= 0.1) %>% 
  mutate(sval = -1*log10(sval)) %>% 
  arrange(Chromosome, Position) %>% 
  inner_join(snps, by = c("SNP", "Chromosome", "Position"))

pA <- ggplot(assoc, aes(x = ps, y = sval, colour = chr_color)) + theme_bw() + 
  geom_vline(aes(xintercept = Bound), bounds, colour = "grey80") + 
  geom_point(size = 1) + 
  geom_hline(yintercept = -log10(0.1), linetype = 2, colour = "red") +
  scale_colour_manual(values = c("grey60", "black")) + 
  scale_x_continuous(breaks = ticks$Tick, labels = ticks$Chromosome) + 
  labs(colour = "", x = "", y = expression(-log[10]*"(s-value)"), tag = "(a)") + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "none", 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 10))


# Phenotypes per SNP ------------------------------------------------------
pB <- get_n_significant_conditions(m2, thresh = 0.1) %>%
  enframe() %>%
  filter(value > 0) %>%
  ggplot(., aes(x = value)) + theme_classic() + 
    geom_histogram(binwidth = 1, colour = "black", fill = "aquamarine", 
                   alpha = 0.8) + 
    labs(x = "# significant phenotypes", y = "# SNPs", tag = "(e)") + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))


# Dominant covariance component -------------------------------------------
# ED1 (relative frequency = 36.2%); ED2 (relative frequency = 11.7%)
pi_mat <- matrix(m2$fitted_g$pi[-1], ncol = length(m2$fitted_g$Ulist), 
                 nrow = length(m2$fitted_g$grid), byrow = TRUE)
colnames(pi_mat) <- names(m2$fitted_g$Ulist)
colSums(pi_mat) %>% sort(decreasing = TRUE) %>% head()

cmat <- cov2cor(m2$fitted_g$Ulist$ED_1)
rownames(cmat) <- colnames(cmat) <- phenotypes
cmat[upper.tri(cmat, diag = TRUE)] <- NA

pC <- as_tibble(cmat, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = phenotypes, ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(phenotypes), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() + 
    geom_tile() +  scale_x_discrete(position = "top") + 
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 90, colour = fills[-1], size = 7), 
          axis.text.y = element_text(colour = rev(fills)[-1], size = 7), 
          legend.position = "inside", 
          legend.position.inside = c(0.275, 0.05), 
          legend.direction = "horizontal", 
          plot.subtitle = element_text(size = 8), 
          legend.title = element_text(size = 7), 
          legend.text = element_text(size = 6), 
          legend.key.width = unit(0.05, "npc"), 
          legend.key.height = unit(0.025, "npc")) + 
    labs(x = "", y = "", fill = expression(r^2), 
         subtitle = "relative frequency = 36.2%", tag = "(b)")


# First eigenvector -------------------------------------------------------
vold <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))$v
u <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))$u
d <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))$d
v <- vold[order(colnames(m2$result$lfsr)), ]

pD <- tibble(Phenotype = phenotypes, 
             V = v[, 1]) %>%
  # mutate(Variable = str_remove(Phenotype, "_X[01]\\.[0-9]{2,3}_X[01]\\.[0-9]{1,3}")) %>%
  ggplot(., aes(x = Phenotype, y = V)) + theme_classic() + 
    geom_col(aes(fill = Phenotype), colour = "black") + 
    labs(x = "", y = "Loadings", fill = "", 
         subtitle = paste0("Eigenvector 1, (PVE = ", 
                           round(d[1]^2/sum(d^2), 2)*100, "%)"), 
         tag = "(c)") + 
    scale_fill_manual(values = fills) +
    guides(fill = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = fills2, size = 7), 
          axis.text.y = element_text(size = 7), 
          axis.title = element_text(size = 8), 
          plot.subtitle = element_text(size = 8))


# Sharing by sign ---------------------------------------------------------
m_data <- read_rds("data/gemma/mash_data.rds")
m_data[-3] <- lapply(m_data[-3], function(m) m[rownames(m2$result$lfsr), ])
pm_mash_beta <- m2$result$PosteriorMean*m_data$s_hat
pm_mash_beta_mag <- pm_mash_beta[rowSums(m2$result$lfsr < 0.1) > 0, ]
lfsr_mash <- m2$result$lfsr[rowSums(m2$result$lfsr < 0.1) > 0, ]
shared_fold_size <- matrix(NA, nrow = ncol(lfsr_mash), ncol = ncol(lfsr_mash))
colnames(shared_fold_size) <- rownames(shared_fold_size) <- phenotypes

for (i in 1:ncol(lfsr_mash)) {
  for (j in 1:ncol(lfsr_mash)) { 
    sig_row <- which(lfsr_mash[, i] < 0.1)
    sig_col <- which(lfsr_mash[, j] < 0.1)
    a <- union(sig_row, sig_col)
    quotient <- pm_mash_beta_mag[a, i]/pm_mash_beta_mag[a, j]
    shared_fold_size[i, j] <- mean(quotient > 0)
  }
}

# shared_fold_size <- shared_fold_size[order(rownames(shared_fold_size)),
#                                      order(colnames(shared_fold_size))]
shared_fold_size[upper.tri(shared_fold_size, diag = TRUE)] <- NA

pE <- as_tibble(shared_fold_size, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = phenotypes, ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(phenotypes), ordered = TRUE)) %>%
  ggplot(aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() + 
    geom_tile() + labs(x = "", y = "", fill = "", tag = "(d)") + 
    scale_fill_distiller(type = "seq", palette = "Greens", 
                         limits = c(0, 1), direction = 1, labels = scales::percent) + 
    scale_x_discrete(position = "top") + 
    theme(axis.text.x = element_text(angle = 90, colour = fills[-1], size = 7), 
          axis.text.y = element_text(colour = rev(fills)[-1], size = 7), 
          legend.position = "inside", 
          legend.direction = "horizontal", 
          legend.position.inside = c(0.275, 0.075), 
          legend.text = element_text(size = 6), 
          legend.key.width = unit(0.05, "npc"), 
          legend.key.height = unit(0.025, "npc"))


# Assemble the figure -----------------------------------------------------
lay <- matrix(c(1, 1, 2, 4, 3, 5), ncol = 2, byrow = TRUE)
gp <- arrangeGrob(ggplotGrob(pA), ggplotGrob(pC), ggplotGrob(pD), 
                  ggplotGrob(pE), ggplotGrob(pB), layout_matrix = lay)
ggsave("figures/final/Figure_5.pdf", gp, width = 180, height = 180, units = "mm", dpi = 600)
