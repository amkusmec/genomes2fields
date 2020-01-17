library(tidyverse)
library(purrrlyr)
library(BGLR)
library(parallel)


# Prepare the data --------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
K <- sommer::A.mat(snps$GD - 1)
D <- sommer::D.mat(snps$GD - 1)

EVA <- eigen(K)
EVD <- eigen(D)

pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(., rownames = "PedigreeNew")

pheno <- read_rds("data/phenotype/glmm_rxn_norm_parameters.rds") %>%
  spread(Parameter, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew)) %>%
  as.data.frame()


# Assign subpopulations ---------------------------------------------------
# Based on fastSTRUCTURE results for K = 8 based on the deltaK method
struct <- read.table("data/structure/g2f_hyb_rep1.8.meanQ", header = FALSE) %>%
  as_tibble() %>%
  mutate(PedigreeNew = rownames(snps$GD)) %>%
  gather(Population, Admixture, V1:V8) %>%
  arrange(desc(Admixture)) %>%
  distinct(PedigreeNew, .keep_all = TRUE)
groups <- struct$Population
names(groups) <- struct$PedigreeNew


# Select stratified training sets -----------------------------------------
sizes <- floor(nrow(pheno)*seq(0.1, 0.9, 0.1))
source("src/st_random.R")
set.seed(687276)

sets <- lapply(sizes, function(i) {
  lapply(1:100, function(j) {
    st_random(rownames(snps$GD), groups, i)
  })
})
actual_size <- sapply(sets, function(l) length(l[[1]]))
names(sets) <- make.names(actual_size)


# Set-up for parallel processing ------------------------------------------
scenarios <- tibble(Phenotype = rep(names(pheno)[2:6], each = length(sizes)*100), 
                    Size = rep(rep(actual_size, each = 100), times = 5), 
                    Replicate = rep(1:100, times = length(sizes)*5))
cl <- makeCluster(detectCores())
clusterEvalQ(cl, { library(tidyverse); library(purrrlyr); library(BGLR) })
clusterExport(cl, list("pheno", "EVA", "EVD", "sets"))
idx <- rep_len(1:detectCores(), length.out = nrow(scenarios)) %>% sort()
scenarios <- split(scenarios, idx)


# Single trait prediction (A model) ---------------------------------------
pred_A <- parLapply(cl, scenarios, function(s) {
  by_row(s, function(r) {
    masked <- pheno[[r$Phenotype[1]]]
    idx <- which(pheno$PedigreeNew %in% sets[[make.names(r$Size[1])]][[r$Replicate[1]]])
    masked[-idx] <- NA
    fm <- BGLR(y = masked, 
               ETA = list(fixed = list(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, 
                                       data = pheno, model = "FIXED"), 
                          A = list(V = EVA$vectors, d = EVA$values, model = "RKHS")), 
               nIter = 12000, burnIn = 2000, 
               saveAt = paste0("data/bglr/trash/A_", r$Phenotype[1], "_", 
                               r$Size[1], "_", r$Replicate[1], "_"))
    yHat <- fm$yHat
    yHat[idx] <- NA
    yHat
  }, .to = "Prediction")
}) %>%
  bind_rows() %>%
  unnest(Prediction) %>%
  mutate(PedigreeNew = rep(pheno$PedigreeNew, times = 5*9*100))
write_rds(pred_A, "data/phenotype/st_A_bglr.rds")

pred_AD <- parLapply(cl, scenarios, function(s) {
  by_row(s, function(r) {
    masked <- pheno[[r$Phenotype[1]]]
    idx <- which(pheno$PedigreeNew %in% sets[[make.names(r$Size[1])]][[r$Replicate[1]]])
    masked[-idx] <- NA
    fm <- BGLR(y = masked, 
               ETA = list(fixed = list(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, 
                                       data = pheno, model = "FIXED"), 
                          A = list(V = EVA$vectors, d = EVA$values, model = "RKHS"), 
                          D = list(V = EVD$vectors, d = EVD$values, model = "RKHS")), 
               nIter = 12000, burnIn = 2000, 
               saveAt = paste0("data/bglr/trash/AD_", r$Phenotype[1], "_", 
                               r$Size[1], "_", r$Replicate[1], "_"))
    yHat <- fm$yHat
    yHat[idx] <- NA
    yHat
  }, .to = "Prediction")
}) %>%
  bind_rows() %>%
  unnest(Prediction) %>%
  mutate(PedigreeNew = rep(pheno$PedigreeNew, times = 5*9*100))
write_rds(pred_AD, "data/phenotype/st_AD_bglr.rds")

stopCluster(cl)
