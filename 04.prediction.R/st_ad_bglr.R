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

pheno <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)
idx <- match(pcs$PedigreeNew, pheno$PedigreeNew)
pheno <- pheno[idx, ]

pheno <- bind_cols(pheno, dplyr::select(pcs, -PedigreeNew)) %>%
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
names(sets) <- make.names(actual_size)


# Set-up for parallel processing ------------------------------------------
scenarios <- tibble(Phenotype = rep(names(pheno)[2:8], each = length(sizes)*100), 
                    Size = rep(rep(actual_size, each = 100), times = 7), 
                    Replicate = rep(1:100, times = length(sizes)*7))
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
               ETA = list(fixed = list(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, 
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
  mutate(PedigreeNew = rep(pheno$PedigreeNew, times = 7*9*100))
write_rds(pred_A, "data/phenotype/st_A_bglr.rds")

pred_AD <- parLapply(cl, scenarios, function(s) {
  by_row(s, function(r) {
    masked <- pheno[[r$Phenotype[1]]]
    idx <- which(pheno$PedigreeNew %in% sets[[make.names(r$Size[1])]][[r$Replicate[1]]])
    masked[-idx] <- NA
    fm <- BGLR(y = masked, 
               ETA = list(fixed = list(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, 
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
  mutate(PedigreeNew = rep(pheno$PedigreeNew, times = 7*9*100))
write_rds(pred_AD, "data/phenotype/st_AD_bglr.rds")

stopCluster(cl)



# Compute correlations ----------------------------------------------------
cors_A <- split(pred_A, list(pred_A$Phenotype, pred_A$Size, pred_A$Replicate)) %>%
  map_df(function(df) {
    y <- pheno[[which(names(pheno) == df$Phenotype[1])]]
    tibble(Phenotype = df$Phenotype[1], 
           Size = df$Size[1], 
           Replicate = df$Replicate[1], 
           R = cor(df$Prediction, y, method = "pearson", use = "complete.obs"))
  }) %>%
  mutate(Model = "A") %>%
  select(Model, Phenotype:R)

cors_AD <- split(pred_AD, list(pred_AD$Phenotype, pred_AD$Size, pred_AD$Replicate)) %>%
  map_df(function(df) {
    y <- pheno[[which(names(pheno) == df$Phenotype[1])]]
    tibble(Phenotype = df$Phenotype[1], 
           Size = df$Size[1], 
           Replicate = df$Replicate[1], 
           R = cor(df$Prediction, y, method = "pearson", use = "complete.obs"))
  }) %>%
  mutate(Model = "A+D") %>%
  select(Model, Phenotype:R)

cors <- bind_rows(cors_A, cors_AD) %>%
  mutate(Model = factor(Model), 
         Size = factor(Size, levels = sort(unique(Size)), ordered = TRUE))

ggplot(cors, aes(x = Size, y = R)) + theme_bw() +
  geom_boxplot(aes(fill = Model), outlier.shape = 1, outlier.colour = "red", 
               alpha = 0.7, notch = TRUE) +
  facet_wrap(~ Phenotype) +
  scale_fill_manual(values = c("A" = "skyblue", "A+D" = "palegreen")) +
  labs(x = "Training Set Size", y = "Prediction Accuracy") +
  theme(panel.grid.major.x = element_blank())
ggsave("figures/prediction/st_ad_accuracy.pdf", width = 9, height = 6, 
       units = "in", dpi = 300)


# Distribution of differences between A and A+D ---------------------------
cors %>%
  spread(Model, R) %>%
  mutate(Difference = `A+D` - A) %>%
  ggplot(., aes(x = Size, y = Difference)) + theme_bw() +
    geom_boxplot(outlier.shape = 1, outlier.colour = "red", notch = TRUE) +
    facet_wrap(~ Phenotype) + theme(panel.grid.major.x = element_blank()) + 
    labs(x = "Training Set Size", y = expression(r[A+D]-r[A]))
ggsave("figures/prediction/st_ad_improvement.pdf", width = 9, height = 6, 
       units = "in", dpi = 300)

cors %>%
  spread(Model, R) %>%
  mutate(Difference = (`A+D` - A)/A) %>%
  filter(abs(Difference) <= 1) %>%
  ggplot(., aes(x = Size, y = Difference)) + theme_bw() +
    geom_boxplot(outlier.shape = 1, outlier.colour = "red", notch = TRUE) +
    facet_wrap(~ Phenotype) + theme(panel.grid.major.x = element_blank()) + 
    labs(x = "Training Set Size", y = "% Improvement") +
    scale_y_continuous(labels = scales::percent)
ggsave("figures/prediction/st_ad_pimprovement.pdf", width = 9, height = 6, 
       units = "in", dpi = 300)
