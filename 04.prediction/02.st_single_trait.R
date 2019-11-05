library(tidyverse)
library(QGenTools)
library(rrBLUP)


# Prepare the data --------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
K <- A.mat(snps$GD)

pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(., rownames = "PedigreeNew")

pheno <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew))
pheno <- as.data.frame(pheno)


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
      # ans <- mixed.solve(y = pheno$Masked, X = cbind(matrix(1, ncol = 1, nrow = nrow(pheno)), 
      #                                                       pheno[, 9:17]), K = K)
      ans <- kin.blup(pheno, "PedigreeNew", "Masked", K = K, 
                      covariate = paste0("PC", 1:9))
      # c(cor(pheno[[p]][-idx], ans$u[-idx], method = "pearson"), 
      #   cor(pheno[[p]][-idx], ans$u[-idx], method = "kendall"))
      u <- rep(NA, nrow(pheno))
      u[-idx] <- ans$pred[-idx]
      u
    })
    # rownames(temp) <- c("Pearson", "Kendall")
    # colnames(temp) <- colnames(pheno)[2:8]
    # as_tibble(t(temp), rownames = "Phenotype") %>%
    #   mutate(Replicate = r)
    colnames(temp) <- colnames(pheno)[2:8]
    temp
  })
  names(rr) <- make.names(1:100)
  rr
})
names(pred) <- make.names(sizes)

write_rds(pred, "data/phenotype/st_single_trait_predictions.rds")

# Apply Fisher's transform to calculate means and standard errors
cors <- map_df(pred, function(l) {
  map_df(l, function(m) {
    tt <- sapply(1:ncol(m), function(x) {
      c(cor(m[, x], pheno[[x + 1]], method = "pearson", use = "complete.obs"), 
        cor(m[, x], pheno[[x + 1]], method = "kendall", use = "complete.obs"))
    }) %>%
      t()
    rownames(tt) <- names(pheno)[2:8]
    colnames(tt) <- c("Pearson", "Kendall")
    as_tibble(tt, rownames = "Phenotype")
  })
}) %>%
  mutate(Size = rep(names(pred), each = 700) %>% str_remove(., "X") %>% as.integer(), 
         Replicate = rep(rep(1:100, each = 7), times = length(sizes))) %>%
  gather(Measure, Value, Pearson:Kendall) %>%
  select(Size, Replicate, Phenotype, Measure, Value)

cors %>%
  mutate(Measure = factor(Measure), 
         Size = factor(Size)) %>%
  ggplot(., aes(x = Size, y = Value)) + theme_bw() + facet_wrap(~ Phenotype) + 
    geom_boxplot(aes(fill = Measure), outlier.colour = "red", outlier.shape = 1) +
    scale_fill_manual(values = c("Pearson" = "palegreen", "Kendall" = "skyblue")) + 
    labs(x = "Training Set Size", y = "Prediction Accuracy")
ggsave("figures/prediction/st_single_trait_accuracy.pdf", width = 8, height = 6, 
       units = "in", dpi = 300)
