library(tidyverse)


# Load the data -----------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
gmap <- snps$GM
snps <- snps$GD
taxa <- rownames(snps)
gc()

pca <- read_rds("data/gbs/pca_covariates.rds")

rxn <- read_rds("data/phenotype/rn_rxn_norm_parameters.rds") %>%
  filter(PedigreeNew %in% taxa) %>%
  spread(Variable, Value)


# Convert the genotypes to BIMBAM format ----------------------------------
bimbam <- gmap %>% 
  select(SNP, Alleles) %>%
  separate(Alleles, c("Allele1", "Allele2"), sep = ",", remove = TRUE) %>%
  cbind(., t(snps))
write_csv(bimbam, "data/gemma/g2f_hybrids.txt", col_names = FALSE)


# Convert the phenotypes --------------------------------------------------
idx <- match(taxa, rxn$PedigreeNew)
rxn <- rxn[idx, ]

# Add a dummy phenotype for relatedness estimation
rxn$Dummy <- rnorm(nrow(rxn))
write_tsv(rxn[, -1], "data/gemma/rxn_norm.txt", col_names = FALSE)


# Convert map information -------------------------------------------------
gmap %>%
  select(SNP, Position, Chromosome) %>%
  write_csv(., "data/gemma/snp_info.txt", col_names = FALSE)


# Convert PC covariates ---------------------------------------------------
# When using covariates, GEMMA requires the explicit specification of an intercept
as_tibble(pca) %>%
  mutate(Intercept = 1) %>%
  select(Intercept, PC1:PC9) %>%
  write_tsv(., "data/gemma/covariates.txt", col_names = FALSE)
