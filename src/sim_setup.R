require(tidyverse)
require(AlphaSimR)
source("src/SimParam_addTraitAD_manual.R")


# Load the pre-computed haplotypes and genetic map ------------------------
haplotypes <- read_rds("data/gbs/synthetic_hybrids_haplotypes.rds")

genMap <- read_rds("data/gbs/sim_map.rds") %>%
  map(function(v) v - min(v))


# Convert SNP positions to AlphaSim map positions -------------------------
snpMap <- read_rds("data/gbs/add_snps.rds")$GM %>%
  split(., .$Chromosome)
effects_A <- read_rds("data/gemma/norm_snp_mash_all.rds")
effects_D <- read_rds("data/gemma/dominance/norm_dom_snp_mash.rds")
loci <- apply(effects_A$result$lfsr, 2, function(x) which(x <= 0.5))
map_pos <- lapply(loci, function(l) {
  lapply(snpMap, function(df) which(df$SNP %in% names(l)))
})


# Load the error covariance structure -------------------------------------
ans <- read_rds("data/phenotype/genetic_corr.rds")
E_init <- ans$sigma$`u:units` %>% cov2cor()


# Load the variance component estimates -----------------------------------
rxn <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)

est <- map_df(2:8, function(p) {
  varE_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_varE.dat"))[-c(1:400)]
  varUA_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  varUD_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_ETA_D_varU.dat"))[-c(1:400)]
  
  tibble(Phenotype = names(rxn)[p], 
         VarA = varUA_D, 
         VarD = varUD_D, 
         VarE = varE_D)
}) %>%
  group_by(Phenotype) %>%
  summarise_all(mean) %>%
  ungroup()


# Set-up the founder population -------------------------------------------
founderPop <- newMapPop(genMap, haplotypes, inbred = FALSE)
SP <- SimParam$new(founderPop)
for (i in 1:length(loci)) {
  SP$addTraitAD_manual(map_pos[[i]], effects_A$result$PosteriorMean[loci[[i]], i, drop = FALSE], 
                       effects_D$result$PosteriorMean[loci[[i]], i, drop = FALSE], 
                       var = est$VarA[i] + est$VarD[i], useVarA = FALSE)
}
SP$setVarE(varE = est$VarE)
SP$setCorE(corE = E_init)


# Randomly sample 100 base populations ------------------------------------
set.seed(847166)
populations <- sapply(1:100, function(i) sample(701, 400, replace = FALSE))


# Utility function for constrained selection indices ----------------------
# Kempthorne and Nordskog (1959)
constrained <- function(wts, G, P, restrict = NULL) {
  if (is.null(restrict)) {
    solve(P) %*% G %*% matrix(wts, ncol = 1)
  } else {
    Gr <- G[restrict, ]
    (diag(nrow(G)) - solve(P) %*% t(Gr) %*% solve(Gr %*% solve(P) %*% t(Gr)) %*% Gr) %*% solve(P) %*% G %*% matrix(wts, ncol = 1)
  }
}
