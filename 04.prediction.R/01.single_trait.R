library(tidyverse)
library(rrBLUP)


# Prepare the data --------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
K <- tcrossprod(scale(snps$GD, center = TRUE, scale = TRUE))/ncol(snps$GD)
write_rds(K, "data/gbs/add_K.rds")

pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(., rownames = "PedigreeNew")

pheno <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew))


# Cross-validation folds --------------------------------------------------
set.seed(827843)

folds <- sapply(1:100, function(i) {
  rep_len(1:10, length.out = nrow(pheno))[sample(1:nrow(pheno), 
                                                 nrow(pheno), replace = FALSE)]
})


# Single-trait prediction with GBLUP --------------------------------------
pred <- lapply(2:8, function(p) {
  cat("\n", names(pheno)[p], "\n")
  temp <- sapply(1:100, function(s) { # Set
    u <- numeric(nrow(pheno))
    for (f in 1:10) { # Iterate over folds
      cat("Set:", s, ", Fold:", f, "\r")
      pheno$Masked <- pheno[[p]]
      idx <- which(folds[, s] == f)
      pheno$Masked[idx] <- NA
      ans <- mixed.solve(y = pheno$Masked, X = cbind(matrix(1, ncol = 1, nrow = nrow(pheno)), 
                                                     pheno[, 9:17]), K = K)
      u[idx] <- ans$u[idx]
    }
    c(cor(pheno[[p]], u, method = "pearson"), 
      cor(pheno[[p]], u, method = "kendall"))
  })
  rownames(temp) <- c("Pearson", "Kendall")
  colnames(temp) <- as.character(1:100)
  as_tibble(t(temp), rownames = "Replicate") %>%
    mutate(Trait = names(pheno)[p])
})

pred <- bind_rows(pred)
write_rds(pred, "data/phenotype/single_trait_predictions.rds")

pred %>%
  gather(Measure, Value, Pearson:Kendall) %>%
  ggplot(., aes(x = Trait, y = Value, fill = Measure)) +
    theme_classic() + geom_boxplot() + facet_wrap(~ Measure, ncol = 1) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(x = "", y = "Correlation", fill = "")
ggsave("figures/prediction/single_trait_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)


# Single-trait prediction with Gaussian kernel ----------------------------
D <- matrix(0, nrow = nrow(pheno), ncol = nrow(pheno))
for (i in 2:nrow(D)) {
  cat(i, "\r")
  for (j in 1:(i - 1)) {
    D[i, j] <- D[j, i] <- sqrt(sum((snps$GD[i, ] - snps$GD[j, ])^2)/(4*ncol(snps$GD)))
  }
}
rownames(D) <- colnames(D) <- rownames(snps$GD)

pheno <- as.data.frame(pheno)
pred2 <- lapply(2:8, function(p) {
  cat("\n", names(pheno)[p], "\n")
  temp <- sapply(1:100, function(s) { # Set
    u <- numeric(nrow(pheno))
    for (f in 1:10) { # Iterate over folds
      cat("Set:", s, ", Fold:", f, "\r")
      pheno$Masked <- pheno[[p]]
      idx <- which(folds[, s] == f)
      pheno$Masked[idx] <- NA
      ans <- kin.blup(pheno, "PedigreeNew", "Masked", GAUSS = TRUE, K = D, 
                      covariate = colnames(pcs)[-1])
      u[idx] <- ans$g[idx]
    }
    c(cor(pheno[[p]], u, method = "pearson"), 
      cor(pheno[[p]], u, method = "kendall"))
  })
  rownames(temp) <- c("Pearson", "Kendall")
  colnames(temp) <- as.character(1:100)
  as_tibble(t(temp), rownames = "Replicate") %>%
    mutate(Trait = names(pheno)[p])
})

pred2 <- bind_rows(pred2)
write_rds(pred2, "data/phenotype/single_trait_predictions_rkhs.rds")

pred2 %>%
  gather(Measure, Value, Pearson:Kendall) %>%
  ggplot(., aes(x = Trait, y = Value, fill = Measure)) +
    theme_classic() + geom_boxplot() + facet_wrap(~ Measure, ncol = 1) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(x = "", y = "Correlation", fill = "")
ggsave("figures/prediction/single_trait_accuracy_rkhs.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
