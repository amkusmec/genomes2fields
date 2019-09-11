library(tidyverse)



# Load the data -----------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
gmap <- snps$GM
snps <- snps$GD
taxa <- rownames(snps)
gc()


# Calculate the dominance relationship matrix -----------------------------
D <- sommer::D.mat(snps - 1)
write_tsv(as_tibble(D), "data/gemma/dominance/dom_mat.txt", col_names = FALSE)


# Convert the genotypes to dominant coding in BIMBAM format ---------------
snps <- apply(snps, 2, function(x) if_else(near(x, 2), 0, if_else(near(x, 1.5), 0.5, x)))
n_alleles <- apply(snps, 2, function(x) length(unique(x)))
bimbam <- gmap %>%
  select(SNP, Alleles) %>%
  separate(Alleles, c("Allele1", "Allele2"), sep = ",", remove = TRUE) %>%
  cbind(., t(snps))
write_csv(bimbam, "data/gemma/dominance/g2f_hybrids_dom.txt", col_names = FALSE)


# Compute PCA for population structure covariates -------------------------
pca <- prcomp(snps, center = TRUE, scale = TRUE)
write_rds(pca, "data/gbs/g2f_hybrids_dom_pca.rds")

as_tibble(pca$x[, 1:7]) %>%
  mutate(Intercept = 1) %>%
  select(Intercept, PC1:PC7) %>%
  write_tsv(., "data/gemma/dominance/covariates_dom.txt", col_names = FALSE)
