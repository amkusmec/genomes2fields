library(tidyverse)


# Keep only hybrids that will be used for GxE analysis --------------------
taxa <- rownames(read_rds("data/gbs/add_snps.rds")$GD)

blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- rep(names(blue), times = sapply(blue, function(x) nrow(x$blue)))
blue <- blue %>%
  map_df(function(x) x$blue) %>%
  mutate(Site = sites) %>%
  filter(PedigreeNew %in% taxa) %>%
  select(Site, PedigreeNew, BLUE)


# Across-site correlations ------------------------------------------------
sites <- unique(sites) %>% sort()
n_sites <- length(sites)
cor_mat <- matrix(0, nrow = n_sites, ncol = n_sites)
diag(cor_mat) <- NA
sig_mat <- cor_mat

for (i in 1:(n_sites - 1)) {
  for (j in (i + 1):n_sites) {
    # Above-diagonal = number of common hybrids
    temp <- blue %>%
      filter(Site == sites[i] | Site == sites[j]) %>%
      spread(Site, BLUE, fill = NA)
    temp <- temp[!is.na(temp[[2]]) & !is.na(temp[[3]]), ]
    cor_mat[i, j] <- nrow(temp)
    
    # Below-diagonal = Kendall's tau b
    cor_mat[j, i] <- cor(temp[[2]], temp[[3]], method = "kendall")
    
    # Significance test
    sig_mat[j, i] <- cor.test(temp[[2]], temp[[3]], method = "kendall")$p.value
  }
}

dimnames(cor_mat) <- dimnames(sig_mat) <- list(sites, sites)


# Plot the correlation matrix ---------------------------------------------
cor_mat_lower <- cor_mat
cor_mat_lower[upper.tri(cor_mat_lower, diag = FALSE)] <- NA
cor_mat_lower <- cor_mat_lower %>%
  as_tibble(rownames = "Site1") %>%
  gather(Site2, Tau, -Site1, na.rm = TRUE)

cor_mat_upper <- cor_mat
cor_mat_upper[lower.tri(cor_mat_upper, diag = FALSE)] <- NA
cor_mat_upper <- cor_mat_upper %>%
  as_tibble(rownames = "Site1") %>%
  gather(Site2, N_hyb, -Site1, na.rm = TRUE)

sig_mat[upper.tri(sig_mat, diag = FALSE)] <- NA
sig_mat <- as_tibble(sig_mat, rownames = "Site1") %>%
  gather(Site2, P, -Site1, na.rm = TRUE) %>%
  mutate(Q = p.adjust(P, method = "fdr"), 
         Sig = if_else(Q <= 0.05, "*", ""))

ggplot() + theme_classic() + 
  geom_tile(aes(x = Site2, y = Site1, fill = Tau), cor_mat_lower) +
  geom_text(aes(x = Site2, y = Site1, label = N_hyb), cor_mat_upper, size = 2) +
  geom_text(aes(x = Site2, y = Site1, label = Sig), sig_mat, size = 4) +
  scale_fill_distiller(type = "div", palette = "RdBu", direction = 1) +
  labs(x = "", y = "", fill = expression(tau[b])) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/single/tau_correlation.pdf", width = 12, height = 12, 
       units = "in", dpi = 300)


cor_mat_upper_sub <- cor_mat_upper %>%
  filter(str_detect(Site1, "NEH"), str_detect(Site2, "NEH"))
cor_mat_lower_sub <- cor_mat_lower %>%
  filter(str_detect(Site1, "NEH"), str_detect(Site2, "NEH"))
sig_mat_sub <- sig_mat %>%
  filter(str_detect(Site1, "NEH"), str_detect(Site2, "NEH"))
ggplot() + theme_classic() + 
  geom_tile(aes(x = Site2, y = Site1, fill = Tau), cor_mat_lower_sub) +
  geom_text(aes(x = Site2, y = Site1, label = N_hyb), cor_mat_upper_sub, size = 5) +
  geom_text(aes(x = Site2, y = Site1, label = Sig), sig_mat_sub, size = 6) + 
  scale_fill_distiller(type = "seq", palette = "Greens", direction = 1) + 
  labs(x = "", y = "", fill = expression(tau[b])) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/single/NEH_tau.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Across-environment correlation (Pearson) --------------------------------
cor_mat2 <- matrix(0, nrow = n_sites, ncol = n_sites)
diag(cor_mat2) <- NA

for (i in 1:(n_sites - 1)) {
  for (j in (i + 1):n_sites) {
    # Above-diagonal = number of common hybrids
    temp <- blue %>%
      filter(Site == sites[i] | Site == sites[j]) %>%
      spread(Site, BLUE, fill = NA)
    temp <- temp[!is.na(temp[[2]]) & !is.na(temp[[3]]), ]
    cor_mat2[i, j] <- nrow(temp)
    
    # Below-diagonal = Pearson's r
    cor_mat2[j, i] <- cor(temp[[2]], temp[[3]], method = "pearson")
  }
}

dimnames(cor_mat2) <- list(sites, sites)


# Plot the correlation matrix (Pearson) -----------------------------------
cor_mat2_lower <- cor_mat2
cor_mat2_lower[upper.tri(cor_mat2_lower, diag = FALSE)] <- NA
cor_mat2_lower <- cor_mat2_lower %>%
  as_tibble(rownames = "Site1") %>%
  gather(Site2, r, -Site1, na.rm = TRUE)

cor_mat2_upper <- cor_mat2
cor_mat2_upper[lower.tri(cor_mat2_upper, diag = FALSE)] <- NA
cor_mat2_upper <- cor_mat2_upper %>%
  as_tibble(rownames = "Site1") %>%
  gather(Site2, N_hyb, -Site1, na.rm = TRUE)

ggplot() + theme_classic() + 
  geom_tile(aes(x = Site2, y = Site1, fill = r), cor_mat2_lower) +
  geom_text(aes(x = Site2, y = Site1, label = N_hyb), cor_mat2_upper, size = 2) +
  scale_fill_distiller(type = "div", palette = "RdBu", direction = 1) +
  labs(x = "", y = "", fill = expression(r)) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/single/pearson_correlation.pdf", width = 12, height = 12, 
       units = "in", dpi = 300)


# Hierarchical clustering on the correlation matrix -----------------------
# Make the correlation matrix symmetric
diag(cor_mat) <- 1
for (i in 1:(ncol(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    cor_mat[i, j] <- cor_mat[j, i]
  }
}

cor_mat[which(is.na(cor_mat))] <- 0
clust <- hclust(as.dist(1 - abs(cor_mat)), method = "ward.D")

## Label by year/state
yr_labs <- clust$labels %>% str_replace(., "[A-Z]{2}H[0-9]_", "") %>%
  factor() %>% as.integer()
yr_cols <- WGCNA::labels2colors(yr_labs)

st_labs <- clust$labels %>% str_replace(., "H[1-4]_201[4-7]", "") %>%
  factor() %>% as.integer()
st_cols <- WGCNA::labels2colors(st_labs)

pdf("figures/single/dendro_tau.pdf", width = 11, height = 6)
WGCNA::plotDendroAndColors(clust, cbind(yr_cols, st_cols), 
                           c("Year", "State"), dendroLabels = NULL, 
                           addGuide = TRUE, main = expression(tau[b]))
dev.off()


# Hierarchical clustering on the correlation matrix (Pearson) -------------
# Make the correlation matrix symmetric
diag(cor_mat2) <- 1
for (i in 1:(ncol(cor_mat2) - 1)) {
  for (j in (i + 1):ncol(cor_mat2)) {
    cor_mat2[i, j] <- cor_mat2[j, i]
  }
}

cor_mat2[which(is.na(cor_mat2))] <- 0
clust2 <- hclust(as.dist(cor_mat2), method = "ward.D")

## Label by year/state
yr_labs2 <- clust2$labels %>% str_replace(., "[A-Z]{2}H[0-9]_", "") %>%
  factor() %>% as.integer()
yr_cols2 <- WGCNA::labels2colors(yr_labs2)

st_labs2 <- clust2$labels %>% str_replace(., "H[1-4]_201[4-7]", "") %>%
  factor() %>% as.integer()
st_cols2 <- WGCNA::labels2colors(st_labs2)

pdf("figures/single/dendro_pearson.pdf", width = 11, height = 6)
WGCNA::plotDendroAndColors(clust2, cbind(yr_cols2, st_cols2), 
                           c("Year", "State"), dendroLabels = NULL, 
                           addGuide = TRUE, main = "r")
dev.off()


library(dendextend)
cor_cophenetic(clust, clust2)
tanglegram(clust, clust2)

fk <- find_k(clust, krange = 2:44)
pam_labs <- WGCNA::labels2colors(fk$pamobject$clustering)
WGCNA::plotDendroAndColors(clust, pam_labs, c(""), dendroLabels = NULL, 
                           addGuide = TRUE, main = expression(tau[b]))
