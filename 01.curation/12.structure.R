library(tidyverse)


# Load the SNP matrix -----------------------------------------------------
hyb <- read_rds("data/gbs/add_snps.rds")
GD <- hyb$GD


# Dissimilarity matrix ----------------------------------------------------
dis <- matrix(0, nrow = nrow(GD), ncol = nrow(GD))
m <- 1/ncol(GD)
for (i in 2:nrow(dis)) {
  for (j in 1:(i - 1)) {
    dis[i, j] <- dis[j, i] <- sqrt(sum(m*(GD[i, ] - GD[j, ])^2))
  }
}
write_rds(dis, "data/gbs/dissim.rds")


# Identify reciprocal hybrids ---------------------------------------------
idx <- which(dis == 0, arr.ind = TRUE)
idx <- idx[idx[,1] != idx[,2], ]
idx <- idx[-c(9, 13, 19:22, 24:34), ]
hybrids <- tibble(H1 = rownames(GD)[idx[, 1]], 
                  H2 = rownames(GD)[idx[, 2]])

# Remove from the genotype matrix
idx_gbs <- which(rownames(GD) %in% hybrids$H2)
GD <- GD[-idx_gbs, ]
maf <- apply(GD, 2, function(x) sum(x)/(2*length(x)))
GD <- GD[, maf >= 0.025 | maf <= 1 - 0.025]
hyb$GM <- hyb$GM[maf >= 0.025 | maf <= 1 - 0.025, ]
hyb$GD17 <- hyb$GD17[, maf >= 0.025 | maf <= 1 - 0.025]
hyb$GD <- GD
write_rds(hyb, "data/gbs/add_snps.rds")

# Change in the yield tables
yield <- read_rds("data/phenotype/yield_agron0.rds")
idx_yield1 <- which(yield$PedigreeNew %in% hybrids$H2)
idx_yield2 <- sapply(idx_yield1, function(x) {
  match(yield$PedigreeNew[x], hybrids$H2)
})
yield$PedigreeNew[idx_yield1] <- hybrids$H1[idx_yield2]
write_rds(yield, "data/phenotype/yield_agron0.rds")

yield <- read_rds("data/phenotype/yield_agron40.rds")
yield$PedigreeNew[idx_yield1] <- hybrids$H1[idx_yield2]
write_rds(yield, "data/phenotype/yield_agron4n0.rds")

# Remove from dissimilarity matrix
dis <- dis[-idx_gbs, -idx_gbs]
write_rds(dis, "data/gbs/dissim.rds")

# Clean-up
rm(hybrids, idx, yield, hyb, idx_gbs, idx_yield1, idx_yield2, maf)
gc()


# nMDS for population structure -------------------------------------------
nmds <- map(2:15, function(d) {
  MASS::isoMDS(dis, k = d)
})
stress <- sapply(nmds, function(x) x$stress)

tibble(Nd = 2:15, Stress = stress) %>%
  ggplot(., aes(x = Nd, y = Stress)) + theme_classic() +
    geom_point(size = 3) + geom_line() +
    labs(x = "Dimensions", y = "Stress") +
    scale_x_continuous(breaks = 2:15)
ggsave("figures/munge/nmds_stress.pdf", width = 6, height = 4, units = "in", dpi = 300)

for (i in seq_along(nmds)) {
  rownames(nmds[[i]]$points) <- rownames(GD)
  colnames(nmds[[i]]$points) <- paste0("Dim", 1:ncol(nmds[[i]]$points))
}
write_rds(nmds, "data/gbs/nmds.rds")


# PCA for population structure --------------------------------------------
pca <- prcomp(GD, center = TRUE, scale = TRUE)
write_rds(pca, "data/gbs/pca.rds")

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
    labs(x = "PC1 (30.4%)", y = "PC2 (4.9%)")
ggsave("figures/munge/pc1_2.pdf", width = 6, height = 4, units = "in", dpi = 300)

write_rds(pca$x[, 1:10], "data/gbs/pca_covariates.rds")