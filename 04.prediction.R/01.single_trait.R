library(tidyverse)
library(rrBLUP)


# Prepare the data --------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
K <- A.mat(snps$GD - 1)

pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(., rownames = "PedigreeNew")

pheno <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew))
pheno <- as.data.frame(pheno)


# Cross-validation folds --------------------------------------------------
set.seed(827843)

folds <- sapply(1:100, function(i) {
  rep_len(1:10, length.out = nrow(pheno))[sample(1:nrow(pheno), 
                                                 nrow(pheno), replace = FALSE)]
})


# Single-trait prediction with GBLUP --------------------------------------
pred <- lapply(2:8, function(p) {
  cat("\n", names(pheno)[p], "\n")
  sapply(1:100, function(s) { # Set
    u <- numeric(nrow(pheno))
    for (f in 1:10) { # Iterate over folds
      cat("Set:", s, ", Fold:", f, "\r")
      pheno$Masked <- pheno[[p]]
      idx <- which(folds[, s] == f)
      pheno$Masked[idx] <- NA
      # ans <- mixed.solve(y = pheno$Masked, X = cbind(matrix(1, ncol = 1, nrow = nrow(pheno)), 
      #                                                pheno[, 9:17]), K = K)
      ans <- kin.blup(pheno, "PedigreeNew", "Masked", K = K, 
                      covariate = paste0("PC", 1:9))
      u[idx] <- ans$pred[idx]
    }
    # c(cor(pheno[[p]], u, method = "pearson"), 
    #   cor(pheno[[p]], u, method = "kendall"))
    u
  })
  # rownames(temp) <- c("Pearson", "Kendall")
  # colnames(temp) <- as.character(1:100)
  # as_tibble(t(temp), rownames = "Replicate") %>%
  #   mutate(Trait = names(pheno)[p])
})

names(pred) <- names(pheno)[2:8]
# pred <- bind_rows(pred)
write_rds(pred, "data/phenotype/single_trait_predictions.rds")





cors <- map_df(1:length(pred), function(m) {
  tt <- apply(pred[[m]], 2, function(x) {
    c(cor(x, pheno[[m + 1]], method = "pearson", use = "complete.obs"), 
      cor(x, pheno[[m + 1]], method = "kendall", use = "complete.obs"))
  }) %>% t()
  rownames(tt) <- 1:100
  colnames(tt) <- c("Pearson", "Kendall")
  as_tibble(tt, rownames = "Replicate") %>%
    mutate(Phenotype = names(pheno)[m + 1])
}) %>%
  mutate(Replicate = rep(1:100, each = 7)) %>%
  gather(Measure, Value, Pearson:Kendall) %>%
  select(Phenotype, Replicate, Measure, Value)

cors %>%
  mutate(Measure = factor(Measure)) %>%
  ggplot(., aes(x = Phenotype, y = Value)) + theme_classic() + 
    geom_boxplot(aes(fill = Measure), outlier.shape = 1, outlier.colour = "red") +
    scale_fill_manual(values = c("Pearson" = "palegreen", "Kendall" = "skyblue")) + 
    labs(x = "", y = "Prediction Accuracy") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/prediction/single_trait_accuracy.pdf", width = 8, height = 6,
       units = "in", dpi = 300)


# Add weather data --------------------------------------------------------
weather <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$"))
ga <- c("NET_X0_X1.275", "NET_X0.65_X1.15", "PPT_X0.525_X0.7", "PPT_X0.4_X1.4", 
        "TMAX_X0.05_X1.375")
weather <- weather[, c("Site", "PedigreeNew", "BLUE", ga)]

idx <- match(weather$PedigreeNew, pheno$PedigreeNew)
pred_large <- lapply(pred, function(x) x[idx, ])
yield_pred <- sapply(1:100, function(i) {
  pred_large$Hybrid[, i] + pred_large$NET_X0.65_X1.15[, i]*weather[["NET_X0.65_X1.15"]] + 
    pred_large$NET_X0_X1.275[, i]*weather[["NET_X0_X1.275"]] +
    pred_large$PPT_X0.4_X1.4[, i]*weather[["PPT_X0.4_X1.4"]] + 
    pred_large$PPT_X0.525_X0.7[, i]*weather[["PPT_X0.525_X0.7"]] +
    pred_large$TMAX_X0.05_X1.375[, i]*weather[["TMAX_X0.05_X1.375"]]
})

yield_cor <- sapply(1:100, function(i) cor(yield_pred[, i], weather$BLUE))

pred %>%
  gather(Measure, Value, Pearson:Kendall) %>%
  ggplot(., aes(x = Trait, y = Value, fill = Measure)) +
    theme_classic() + geom_boxplot() + facet_wrap(~ Measure, ncol = 1) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(x = "", y = "Correlation", fill = "")
ggsave("figures/prediction/single_trait_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
