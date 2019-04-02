library(tidyverse)


# Load the data -----------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
gmap <- snps$GM
snps <- snps$GD
taxa <- rownames(snps)

pca <- read_rds("data/gbs/pca_covariates.rds")

blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- rep(names(blue), times = sapply(blue, function(x) nrow(x$blue)))
blue <- blue %>%
  map_df(function(x) x$blue) %>%
  mutate(Site = sites) %>%
  filter(PedigreeNew %in% taxa) %>%
  filter(!str_detect(Site, "2017")) %>%
  select(Site, PedigreeNew, BLUE) %>%
  spread(Site, BLUE, fill = NA)


# Convert the genotypes to BIMBAM format ----------------------------------
bimbam <- gmap %>%
  select(SNP, Alleles) %>%
  separate(Alleles, c("Allele1", "Allele2"), sep = ",", remove = TRUE) %>%
  cbind(., t(snps))
write_csv(bimbam, "data/gemma/g2f_hybrids.txt", col_names = FALSE)


# Convert the phenotypes --------------------------------------------------
idx <- match(taxa, blue$PedigreeNew)
blue <- blue[idx, ]

# Add a dummy phenotype for relatedness estimation
blue$Dummy <- rnorm(nrow(blue))
write_tsv(blue[, -1], "data/gemma/yield_blue.txt", col_names = FALSE)


# Convert map information -------------------------------------------------
gmap %>%
  select(SNP, Position, Chromosome) %>%
  write_csv(., "data/gemma/snp_info.txt", col_names = FALSE)


# Convert PC covariates ---------------------------------------------------
# When using covariates, GEMMA requires the explicit specification of an
# intercept.
as_tibble(pca) %>%
  mutate(Intercept = 1) %>%
  select(Intercept, PC1:PC9) %>%
  write_tsv(., "data/gemma/covariates.txt", col_names = FALSE)
