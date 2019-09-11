library(tidyverse)
library(sommer)


# Prepare the data --------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
K <- A.mat(snps$GD - 1)
D <- D.mat(snps$GD - 1)

pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(., rownames = "PedigreeNew")

pheno <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew)) %>%
  mutate(PedigreeNewD = factor(PedigreeNew), 
         PedigreeNew = factor(PedigreeNew)) %>%
  as.data.frame()


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
actual_size <- sapply(sets, function(l) length(l[[1]]))


# Single-trait prediction with GBLUP --------------------------------------
pred_A <- lapply(1:length(sets), function(s) { # Set size
  rr <- lapply(1:100, function(r) { # Replicate
    cat("Size:", sizes[s], ", Replicate:", r, "\r")
    temp <- sapply(2:8, function(p) { # Phenotype
      pheno$Masked <- pheno[[p]]
      idx <- which(pheno$PedigreeNew %in% sets[[s]][[r]])
      pheno$Masked[-idx] <- NA
      ans <- mmer(Masked ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, 
                  random = ~ vs(PedigreeNew, Gu = K), rcov = ~ units, data = pheno, 
                  verbose = FALSE)
      u <- rep(NA, nrow(pheno))
      u[-idx] <- ans$U$`u:PedigreeNew`$Masked[-idx]
      u
    })
    colnames(temp) <- colnames(pheno)[2:8]
    temp
  })
  names(rr) <- make.names(1:100)
  rr
})
names(pred_A) <- make.names(actual_size)

write_rds(pred_A, "data/phenotype/st_single_trait_predictions.rds")


pred_AD <- lapply(1:length(sets), function(s) { # Set size
  rr <- lapply(1:100, function(r) { # Replicate
    cat("Size:", sizes[s], ", Replicate:", r)
    temp <- sapply(2:8, function(p) { # Phenotype
      pheno$Masked <- pheno[[p]]
      idx <- which(pheno$PedigreeNew %in% sets[[s]][[r]])
      pheno$Masked[-idx] <- NA
      ans <- mmer(Masked ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, 
                  random = ~ vs(PedigreeNew, Gu = K) + vs(PedigreeNewD, Gu = D), 
                  rcov = ~ units, data = pheno, verbose = FALSE, tolparinv = 1e-2)
      u <- rep(NA, nrow(pheno))
      if (length(ans) > 0) {
        u[-idx] <- ans$U$`u:PedigreeNew`$Masked[-idx] + 
          ans$U$`u:PedigreeNewD`$Masked[-idx]
      }
      u
    })
    colnames(temp) <- colnames(pheno)[2:8]
    temp
  })
  names(rr) <- make.names(1:100)
  rr
})
names(pred_AD) <- make.names(actual_size)
write_rds(pred_AD, "data/phenotype/st_single_trait_predictions_AD.rds")


# cors <- map_df(pred, function(l) {
#   map_df(l, function(m) {
#     tt <- sapply(1:ncol(m), function(x) {
#       c(cor(m[, x], pheno[[x + 1]], method = "pearson", use = "complete.obs"), 
#         cor(m[, x], pheno[[x + 1]], method = "kendall", use = "complete.obs"))
#     }) %>%
#       t()
#     rownames(tt) <- names(pheno)[2:8]
#     colnames(tt) <- c("Pearson", "Kendall")
#     as_tibble(tt, rownames = "Phenotype")
#   })
# }) %>%
#   mutate(Size = rep(names(pred), each = 700) %>% str_remove(., "X") %>% as.integer(), 
#          Replicate = rep(rep(1:100, each = 7), times = length(sizes))) %>%
#   gather(Measure, Value, Pearson:Kendall) %>%
#   select(Size, Replicate, Phenotype, Measure, Value)
# 
# cors %>%
#   mutate(Measure = factor(Measure), 
#          Size = factor(Size)) %>%
#   ggplot(., aes(x = Size, y = Value)) + theme_bw() + facet_wrap(~ Phenotype) + 
#   geom_boxplot(aes(fill = Measure), outlier.colour = "red", outlier.shape = 1) +
#   scale_fill_manual(values = c("Pearson" = "palegreen", "Kendall" = "skyblue")) + 
#   labs(x = "Training Set Size", y = "Prediction Accuracy")
# ggsave("figures/prediction/st_single_trait_accuracy.pdf", width = 8, height = 6, 
#        units = "in", dpi = 300)

cors_A <- map_df(pred_A, function(l) {
  map_df(l, function(m) {
    tt <- sapply(1:ncol(m), function(x) {
      cor(m[, x], pheno[[x + 1]], method = "pearson", use = "complete.obs")
    })
    names(tt) <- names(pheno)[2:8]
    enframe(tt, name = "Phenotype")
  })
}) %>%
  mutate(Size = rep(names(pred_A), each = 700) %>% str_remove(., "X") %>% as.integer(), 
         Replicate = rep(rep(1:100, each = 7), times = length(sizes))) %>%
  rename(Value = value) %>%
  dplyr::select(Size, Replicate, Phenotype, Value)

ggplot(cors_A, aes(x = Size, y = Value, group = Size)) + theme_bw() +
  geom_boxplot() + facet_wrap(~ Phenotype)
