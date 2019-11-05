library(tidyverse)
library(parallel)
# library(lme4)


# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  separate(Site, c("Environment", "Year"), sep = "_", remove = FALSE) %>%
  select(BLUE, PedigreeNew, Site, Environment, Year, everything())
net <- names(data)[which(str_detect(names(data), "NET"))]
data <- mutate_at(data, net, function(x) -1*x)
ped_site <- data %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# Response
y <- data$BLUE
y <- y - mean(y)

# Predictors
ga <- read_csv("data/weather/ga_windows2.csv")
vars <- paste(ga$Category, make.names(ga$Start), make.names(ga$End), sep = "_") %>%
  unique()
X <- cbind(model.matrix(~ 0 + Environment + Year, data = data), 
           as.matrix(data[, vars]))
X[, -c(1:30)] <- scale(X[, -c(1:30)], center = TRUE, scale = FALSE)

# Variance-covariance matrix of the BLUEs
source("src/direct_sum.R")
d <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
d <- d[!str_detect(names(d), "2017$")]
d <- mapply(function(x, y) {
  idx <- which(rownames(x$vcov) %in% y$PedigreeNew)
  x$vcov[idx, idx]
}, x = d, y = ped_site)
R <- direct_sum(d)
wts <- 1/diag(R)


# Evaluate all possible models --------------------------------------------
models <- lapply(1:5, function(i) {
  temp <- combn(vars, i)
  if (nrow(temp) < 5) {
    temp <- rbind(temp, matrix(NA, ncol = ncol(temp), nrow = 5 - i))
  }
  return(temp)
})
models <- do.call("cbind", models)

cl <- makeCluster(detectCores())
clusterExport(cl, list("y", "X", "wts"))
bic <- parApply(cl, models, 2, function(v) {
  v <- v[!is.na(v)]
  temp <- X[, c(1:30, which(colnames(X) %in% v))]
  BIC(lm(y ~ 0 + temp, weights = wts)) + 2*1*log(choose(ncol(X) - 30, length(v)))
})
stopCluster(cl)

# Add the model without environmental variables
# N.B.: The ith column of `models` corresponds to the (i+1)th element of `bic`.
bic <- c(BIC(lm(y ~ 0 + X[, 1:30], weights = wts)) + 
           2*1*log(choose(ncol(X) - 30, 0)), bic)

selected <- models[, which.min(bic) - 1]


# Post-processing before reaction norm calculation ------------------------
data <- data[, c("Site", "Environment", "Year", "PedigreeNew", "BLUE", selected)]

# Save for future use
write_rds(data, "data/phenotype/yield_blue_final2.rds")

# A site-wise table for post hoc analysis of single-site results
data %>%
  group_by(Site) %>%
  summarise_at(c("BLUE", selected), median) %>%
  ungroup() %>%
  write_rds(., "data/phenotype/yield_site_final2.rds")


# # Compute the joint reaction norms ----------------------------------------
# form <- as.formula(paste0("BLUE ~ 1 + Year + Environment + ", 
#                           paste(selected, collapse = " + "), 
#                           " + (1 + ", 
#                           paste(selected, collapse = " + "), 
#                           " | PedigreeNew)"))
# data[, 6:10] <- scale(data[, 6:10], center = TRUE, scale = FALSE)
# rxn <- lmer(form, data = data, REML = TRUE, weights = wts, 
#             control = lmerControl(optCtrl = list(maxfun = 1e6)))
# 
# eff <- ranef(rxn)[[1]] %>%
#   as_tibble(rownames = "PedigreeNew") %>%
#   rename(Hybrid = `(Intercept)`)
# rr <- data %>%
#   select(PedigreeNew) %>%
#   mutate(E = residuals(rxn)) %>%
#   group_by(PedigreeNew) %>%
#   summarise(VarE = log(var(E))) %>%
#   ungroup()
# eff <- inner_join(eff, rr, by = "PedigreeNew")
# 
# eff %>%
#   gather(Phenotype, Value, -PedigreeNew) %>%
#   ggplot(., aes(x = Value)) + theme_classic() +
#     geom_density() + facet_wrap(~ Phenotype, scales = "free")
# 
# write_rds(eff, "data/phenotype/rxn_norm_parameters.rds")
