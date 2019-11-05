library(tidyverse)
library(BGLR)

rxn <- read_rds("data/phenotype/rn_rxn_norm_parameters.rds") %>%
  spread(Variable, Value)

snps <- read_rds("data/gbs/add_snps.rds")
K <- sommer::A.mat(snps$GD - 1)
D <- sommer::D.mat(snps$GD - 1)

idx <- match(rownames(K), rxn$PedigreeNew)
rxn <- rxn[idx, ]
pca <- read_rds("data/gbs/pca_covariates.rds")
rxn <- rxn %>%
  inner_join(., as_tibble(pca, rownames = "PedigreeNew"), by = "PedigreeNew")

EVA <- eigen(K)
EVD <- eigen(D)


# Run the models ----------------------------------------------------------
for (p in 2:8) {
  fm_A <- BGLR(y = rxn[[p]], 
               ETA = list(fixed = list(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, 
                                       data = rxn, model = "FIXED"), 
                          A = list(V = EVA$vectors, d = EVA$values, model = "RKHS")), 
               nIter = 12000, burnIn = 2000, 
               saveAt = paste0("data/bglr/h2/eigA_", names(rxn)[p], "_"))
  fm_D <- BGLR(y = rxn[[p]], 
               ETA = list(fixed = list(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, 
                                       data = rxn, model = "FIXED"), 
                          A = list(V = EVA$vectors, d = EVA$values, model = "RKHS"), 
                          D = list(V = EVD$vectors, d = EVD$values, model = "RKHS")), 
               nIter = 12000, burnIn = 2000, 
               saveAt = paste0("data/bglr/h2/eigD_", names(rxn)[p], "_"))
}


# Collect the estimates ---------------------------------------------------
phenotypes <- c("Hybrid", "Residual variance", "Drought (early)", 
                "Drought (anthesis)", "Max. temp. (anthesis)", 
                "Min. temp. (season)", "Min. temp. (early)")
est <- map_df(2:8, function(p) {
  varE_A <- scan(paste0("data/bglr/h2/eigA_", names(rxn)[p], "_varE.dat"))[-c(1:400)]
  varU_A <- scan(paste0("data/bglr/h2/eigA_", names(rxn)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  
  varE_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_varE.dat"))[-c(1:400)]
  varUA_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  varUD_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_ETA_D_varU.dat"))[-c(1:400)]
  
  tibble(Phenotype = phenotypes[p - 1], 
         Model = rep(c("A", "A+D"), each = length(varE_A)), 
         h2 = c(varU_A/(varU_A + varE_A), 
                (varUA_D + varUD_D)/(varUA_D + varUD_D + varE_D)))
})

p <- est %>% 
  mutate(Model = factor(Model), 
         Phenotype = factor(Phenotype, levels = unique(Phenotype), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype, y = h2)) + theme_bw() +
    geom_boxplot(aes(fill = Model), outlier.shape = 1, outlier.colour = "red", 
                 position = position_dodge(0.9), alpha = 0.7, notch = TRUE) +
    scale_fill_manual(values = c("A" = "skyblue", "A+D" = "palegreen")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) + 
    labs(x = "", y = expression(h[SNP]^2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.x = element_blank())
ggsave("figures/select/bglr_h2.pdf", plot = p, width = 8, height = 5, 
       units = "in", dpi = 300)
