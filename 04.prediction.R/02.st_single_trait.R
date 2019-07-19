library(tidyverse)
library(QGenTools)
library(rrBLUP)


# Prepare the data --------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
K <- realized_ab(snps$GD)

pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(., rownames = "PedigreeNew")

pheno <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew))


# Assign subpopulations ---------------------------------------------------
# Based on fastSTRUCTURE results for K = 5 based on the deltaK method
struct <- read.table("data/structure/g2f_hyb_rep1.5.meanQ", header = FALSE) %>%
  as_tibble() %>%
  mutate(PedigreeNew = rownames(snps$GD)) %>%
  gather(Population, Admixture, V1:V5) %>%
  arrange(desc(Admixture)) %>%
  distinct(PedigreeNew, .keep_all = TRUE)
struct$Population[struct$PedigreeNew == "Z022E0108/LH82"] <- "V2"
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


# Single-trait prediction with GBLUP --------------------------------------
pred <- lapply(1:length(sets), function(s) { # Set size
  rr <- lapply(1:100, function(r) { # Replicate
    temp <- sapply(2:8, function(p) { # Phenotype
      cat("Size:", sizes[s], ", Replicate:", r, "\r")
      pheno$Masked <- pheno[[p]]
      idx <- which(pheno$PedigreeNew %in% sets[[s]][[r]])
      pheno$Masked[-idx] <- NA
      ans <- mixed.solve(y = pheno$Masked, X = cbind(matrix(1, ncol = 1, nrow = nrow(pheno)), 
                                                            pheno[, 9:17]), K = K)
      c(cor(pheno[[p]][-idx], ans$u[-idx], method = "pearson"), 
        cor(pheno[[p]][-idx], ans$u[-idx], method = "kendall"))
    })
    rownames(temp) <- c("Pearson", "Kendall")
    colnames(temp) <- colnames(pheno)[2:8]
    as_tibble(t(temp), rownames = "Phenotype") %>%
      mutate(Replicate = r)
  }) %>%
    bind_rows() %>%
    mutate(Size = sizes[s])
}) %>%
  bind_rows()

write_rds(pred, "data/phenotype/st_single_trait_predictions.rds")

# Apply Fisher's transform to calculate means and standard errors

pred %>%
  group_by(Phenotype, Size) %>%
  summarise(Mean = mean(0.5*log((1 + Pearson)/(1 - Pearson))), 
            SE = sd(0.5*log((1 + Pearson)/(1 - Pearson)))/sqrt(n())) %>%
  ungroup() %>%
  ggplot(., aes(x = Size, y = Mean, group = Phenotype)) + theme_classic() +
    geom_line(aes(colour = Phenotype)) + 
    geom_pointrange(aes(colour = Phenotype, ymin = Mean - SE, ymax = Mean + SE))
