library(tidyverse)
library(sommer)
library(parallel)


pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(rownames = "PedigreeNew")
rxn <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value) %>%
  inner_join(., pcs, by = "PedigreeNew")

snps <- read_rds("data/gbs/add_snps.rds")
idx <- match(rxn$PedigreeNew, rownames(snps$GD))
snps$GD <- snps$GD[idx, ]
K <- A.mat(snps$GD - 1)

fix1 <- paste0("cbind(", paste(names(rxn)[2:8], collapse = ", "), ")")
fix2 <- paste(paste0("PC", 1:9), collapse = " + ")
fix <- paste0(fix1, " ~ 1 + ", fix2) %>% as.formula()

set.seed(757561)
straps <- sapply(1:1000, function(i) sample(nrow(rxn), nrow(rxn), replace = TRUE))

cl <- makeCluster(40)
clusterEvalQ(cl, library(sommer))
clusterExport(cl, list("rxn", "K", "fix", "straps"))

m <- parLapply(cl, 1:ncol(straps), function(i) {
  rxn2 <- rxn[straps[, i], ]
  K2 <- K[straps[, i], straps[, i]]
  ans <- mmer(fix, random = ~ vs(PedigreeNew, Gu = K2, Gtc = unsm(7)), 
              rcov = ~ vs(units, Gtc = unsm(7)), 
              data = rxn2, date.warning = FALSE)
  list(P = cov(rxn2[, 2:8]), G = ans$sigma$`u:PedigreeNew`, R = ans$sigma$`u:units`)
})

write_rds(m, "data/phenotype/bootstrap_samples.rds")
stopCluster(cl)

P <- lapply(m, function(x) cov2cor(x$P)) %>% abind::abind(along = 3)
G <- lapply(m, function(x) cov2cor(x$G)) %>% abind::abind(along = 3)
R <- lapply(m, function(x) cov2cor(x$R)) %>% abind::abind(along = 3)

temp <- P[, , 1]
temp[upper.tri(temp, diag = TRUE)] <- NA
bounds <- as_tibble(temp, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  dplyr::select(-R) %>%
  purrrlyr::by_row(function(r) {
    tibble(P_lower = quantile(P[r$Phenotype1[1], r$Phenotype2[1], ], probs = 0.025), 
           P_upper = quantile(P[r$Phenotype1[1], r$Phenotype2[1], ], probs = 0.975), 
           G_lower = quantile(G[r$Phenotype1[1], r$Phenotype2[1], ], probs = 0.025), 
           G_upper = quantile(G[r$Phenotype1[1], r$Phenotype2[1], ], probs = 0.975), 
           R_lower = quantile(R[r$Phenotype1[1], r$Phenotype2[1], ], probs = 0.025), 
           R_upper = quantile(R[r$Phenotype1[1], r$Phenotype2[1], ], probs = 0.975))
  }) %>%
  unnest() %>%
  mutate(P = paste0("(", round(P_lower, 2), ",", round(P_upper, 2), ")"), 
         G = paste0("(", round(G_lower, 2), ",", round(G_upper, 2), ")"), 
         R = paste0("(", round(R_lower, 2), ",", round(R_upper, 2), ")"))

write_rds(bounds, "data/phenotype/bootstrap_intervals.rds")
