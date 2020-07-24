### Principal components analysis on the hybrid genotype matrix to estimate 
### population structure.

library(tidyverse)


# Load the SNP matrix -----------------------------------------------------
hyb <- read_rds("data/gbs/add_snps.rds")
GD <- hyb$GD


# Dissimilarity matrix ----------------------------------------------------
### Can be run in the background with data/gbs/dissim.R
# dis <- matrix(0, nrow = nrow(GD), ncol = nrow(GD))
# m <- 1/ncol(GD)
# for (i in 2:nrow(dis)) {
#   for (j in 1:(i - 1)) {
#     dis[i, j] <- dis[j, i] <- sqrt(sum(m*(GD[i, ] - GD[j, ])^2))
#   }
# }
# write_rds(dis, "data/gbs/dissim.rds")
# dis <- read_rds("data/gbs/dissim.rds")
# 
# 
# # nMDS for population structure -------------------------------------------
# nmds <- map(2:15, function(d) {
#   MASS::isoMDS(dis, k = d)
# })
# stress <- sapply(nmds, function(x) x$stress)
# 
# tibble(Nd = 2:15, Stress = stress) %>%
#   ggplot(., aes(x = Nd, y = Stress)) + theme_classic() +
#     geom_point(size = 3) + geom_line() +
#     labs(x = "Dimensions", y = "Stress") +
#     scale_x_continuous(breaks = 2:15)
# ggsave("figures/munge/nmds_stress.pdf", width = 6, height = 4, units = "in", dpi = 300)
# 
# for (i in seq_along(nmds)) {
#   rownames(nmds[[i]]$points) <- rownames(GD)
#   colnames(nmds[[i]]$points) <- paste0("Dim", 1:ncol(nmds[[i]]$points))
# }
# write_rds(nmds, "data/gbs/nmds.rds")


# PCA for population structure --------------------------------------------
### Can be run in the background with data/gbs/pca.R
pca <- prcomp(GD, center = TRUE, scale = TRUE)
write_rds(pca, "data/gbs/pca.rds")
# pca <- read_rds("data/gbs/pca.rds")

var_exp <- pca$sdev^2
var_exp <- var_exp/sum(var_exp)
tibble(PCs = seq_along(var_exp), 
       VExp = var_exp) %>%
  filter(PCs <= 20) %>%
  ggplot(., aes(x = PCs, y = VExp)) + theme_classic() +
    geom_point(size = 3) + geom_line(linetype = 2) +
    geom_hline(yintercept = 0.02, linetype = 3) +
    labs(x = "PC", y = "Variance Explained") +
    scale_y_continuous(labels = scales::percent)
ggsave("figures/munge/scree_plot.pdf", width = 6, height = 4, units = "in", dpi = 300)

as_tibble(pca$x[, 1:2]) %>%
  ggplot(., aes(x = PC1, y = PC2)) + theme_classic() +
    geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "PC1 (17.3%)", y = "PC2 (5.6%)")
ggsave("figures/munge/pc1_2.pdf", width = 6, height = 4, units = "in", dpi = 300)

# All PCs that explain >= 2% of the variance
write_rds(pca$x[, which(var_exp >= 0.02)], "data/gbs/pca_covariates.rds")
