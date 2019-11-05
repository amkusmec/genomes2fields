library(tidyverse)
library(BGLR)
library(parallel)
library(coda)


# Load the data -----------------------------------------------------------
# Phenotype and weather data
data <- read_rds("data/phenotype/yield_blue_final2.rds") %>%
  separate(Site, c("Location", "Year"), sep = "_", remove = FALSE)

ga <- c("NET_X0.025_X0.45", "NET_X0.65_X1.275", "TMAX_X0.775_X0.875", 
        "TMIN_X0.05_X1.5", "TMIN_X0.15_X0.65")


# Center and scale the predictors -----------------------------------------
means <- apply(data[, 7:11], 2, mean)
# sds <- apply(data[, 7:11], 2, sd)
# 
for (i in 7:11) {
#   data[[i]] <- (data[[i]] - means[i - 6])/sds[i - 6]
  data[[i]] <- data[[i]] - means[i - 6]
}


# Compute joint reaction norms --------------------------------------------
# Construct a model object
ETA <- list(fixed = list(~ factor(Site) + NET_X0.025_X0.45 + NET_X0.65_X1.275 + 
                           TMAX_X0.775_X0.875 + TMIN_X0.05_X1.5 + TMIN_X0.15_X0.65, 
                         data = data, model = "FIXED", saveEffects = TRUE), 
            hybrid = list(~ 0 + factor(PedigreeNew), data = data, model = "BRR", 
                          saveEffects = TRUE), 
            NET_X0.025_X0.45 = list(~ 0 + factor(PedigreeNew):NET_X0.025_X0.45, 
                                    data = data, model = "BRR", saveEffects = TRUE), 
            NET_X0.65_X1.275 = list(~ 0 + factor(PedigreeNew):NET_X0.65_X1.275, 
                                    data = data, model = "BRR", saveEffects = TRUE), 
            TMAX_X0.775_X0.875 = list(~ 0 + factor(PedigreeNew):TMAX_X0.775_X0.875, 
                                      data = data, model = "BRR", saveEffects = TRUE), 
            TMIN_X0.05_X1.5 = list(~ 0 + factor(PedigreeNew):TMIN_X0.05_X1.5, 
                                   data = data, model = "BRR", saveEffects = TRUE), 
            TMIN_X0.15_X0.65 = list(~ 0 + factor(PedigreeNew):TMIN_X0.15_X0.65, 
                                    data = data, model = "BRR", saveEffects = TRUE))

# Construct groups for heterogeneous error variances
groups <- as.integer(factor(data$Site))

# Run 4 chains
cl <- makeCluster(4)
clusterEvalQ(cl, library(BGLR))
clusterExport(cl, list("ETA", "groups", "data"))

rxn <- parLapply(cl, 1:4, function(i) {
  BGLR(y = data$BLUE, ETA = ETA, burnIn = 1e5, nIter = 1.1e6, thin = 1e3, 
       saveAt = paste0("data/bglr/reversed_net/rxn_chain", i, "_"), groups = groups)
})

stopCluster(cl)
write_rds(rxn, "data/phenotype/reversed_net_joint_rxn_norm_model.rds")


# Check chain convergence -------------------------------------------------
site_levels <- sort(unique(data$Site))
files <- list.files("data/bglr/reversed_net", "*\\.", full.names = TRUE) %>%
  split(., rep(1:4, each = 15))
mcmc_res_list <- lapply(files, function(f) {
  temp <- lapply(f, function(f2) {
    if (str_detect(f2, "dat$") & str_detect(f2, "fixed")) {
      temp <- read_delim(f2, delim = " ")
    } else if (str_detect(f2, "dat$") & !str_detect(f2, "fixed")) {
      temp <- read_delim(f2, delim = " ", col_names = FALSE)
    } else {
      temp <- readBinMat(f2)
    }

    if (!str_detect(f2, "fixed")) {
      if (str_detect(f2, "varB")) {
        colnames(temp) <- str_replace(f2, "data/bglr/reversed_net/rxn_chain[1-4]_ETA_", "") %>%
          str_replace(., "\\.dat", "")
      } else if (str_detect(f2, "mu")) {
        colnames(temp) <- "mu"
      } else if (str_detect(f2, "varE")) {
        colnames(temp) <- site_levels[as.integer(names(rxn[[1]]$varE))]
      } else {
        comp <- str_replace(f2, "data/bglr/reversed_net/rxn_chain[1-4]_ETA_", "") %>%
          str_replace(., "_b\\.bin", "")
        colnames(temp) <- rxn[[1]]$ETA[[comp]]$colNames %>%
          str_replace(., "factor\\(PedigreeNew\\)", "")
      }
    }

    if (nrow(temp) > 1000) temp <- temp[101:1100, ]

    return(as.matrix(temp))
  })

  mcmc(do.call("cbind", temp), start = 1, end = 1000)
}) %>%
  as.mcmc.list()

# 4,309 parameters total
gdiag <- geweke.diag(mcmc_res_list)
idx <- lapply(gdiag, function(g) which(abs(g$z) < 1.96))
length(reduce(idx, intersect)) # 3,314/4,309 (77%) parameters converged in all chains
length(reduce(idx, union)) # 100% parameters converged in at least one chain
unlist(idx) %>% enframe() %>% count(value) %>% count(n)


# Collect the Gibbs samples -----------------------------------------------
gibbs <- lapply(mcmc_res_list, as.matrix)
gibbs <- do.call("rbind", gibbs)
params <- apply(gibbs, 2, mean) %>%
  enframe()

mus <- filter(params, name %in% ga)
ped <- filter(params, str_detect(name, "/")) %>%
  separate("name", c("PedigreeNew", "Variable"), sep = ":", remove = TRUE) %>%
  mutate(Variable = if_else(is.na(Variable), "Hybrid", Variable)) %>%
  full_join(., mus, by = c("Variable" = "name")) %>%
  mutate(Value = if_else(!is.na(value.y), value.x + value.y, value.x)) %>%
  select(-value.x, -value.y)

# Compute MSE for each hybrid
mu <- mean(gibbs[, "mu"])
model_mat <- cbind(model.matrix(~ 0 + factor(Site) + NET_X0.025_X0.45 + NET_X0.65_X1.275 +
                                  TMAX_X0.775_X0.875 + TMIN_X0.05_X1.5 + TMIN_X0.15_X0.65,
                                data = data)[, -1],
                   model.matrix(~ 0 + factor(PedigreeNew), data = data),
                   model.matrix(~ 0 + factor(PedigreeNew):NET_X0.025_X0.45, data = data),
                   model.matrix(~ 0 + factor(PedigreeNew):NET_X0.65_X1.275, data = data),
                   model.matrix(~ 0 + factor(PedigreeNew):TMAX_X0.775_X0.875, data = data),
                   model.matrix(~ 0 + factor(PedigreeNew):TMIN_X0.05_X1.5, data = data),
                   model.matrix(~ 0 + factor(PedigreeNew):TMIN_X0.15_X0.65, data = data))
colnames(model_mat) <- if_else(str_detect(colnames(model_mat), "factor"),
                               str_replace(colnames(model_mat), "factor\\(PedigreeNew\\)", ""), # %>%
                               #   str_split(., ":") %>%
                               #   sapply(., function(x) paste(rev(x), collapse = ":")),
                               colnames(model_mat))
beta_hat <- params %>%
  slice(1:4261) %>%
  filter(!str_detect(name, "var"))
idx <- match(colnames(model_mat), beta_hat$name)
beta_hat <- beta_hat[idx, ]

y_hat <- drop(mu + model_mat %*% matrix(beta_hat$value, ncol = 1))
ped <- data %>%
  mutate(E = BLUE - y_hat) %>%
  select(PedigreeNew, E) %>%
  group_by(PedigreeNew) %>%
  summarise(Value = log(mean(E^2))) %>%
  ungroup() %>%
  mutate(Variable = "lnMSE") %>%
  select(PedigreeNew, Variable, Value) %>%
  bind_rows(., ped)
write_rds(ped, "data/phenotype/rn_rxn_norm_parameters.rds")

ggplot(ped, aes(x = Value)) + theme_classic() +
  geom_density() + facet_wrap(~ Variable, scales = "free") +
  labs(x = "", y = "Density")

rxn <- ped %>%
  spread(Variable, Value)
rxn_cor <- cor(rxn[, -1])
rxn_cor[upper.tri(rxn_cor, diag = TRUE)] <- NA
sig <- tibble(Phenotype2 = rep(rownames(rxn_cor)[-7], times = seq(6, 1, -1)),
              Phenotype1 = map(1:6, function(i) colnames(rxn_cor)[-c(1:i)]) %>%
                unlist()) %>%
  purrrlyr::by_row(function(r) {
    cor.test(rxn[[r$Phenotype1[1]]], rxn[[r$Phenotype2[1]]])$p.value
  }, .to = "p_val") %>%
  mutate(p_val = p.adjust(p_val, method = "fdr"),
         p_val = if_else(p_val <= 0.05, "*", as.character(NA)))
as_tibble(rxn_cor, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE),
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2)) + theme_classic() +
    geom_tile(aes(fill = R)) + labs(x = "", y = "", fill = "") +
    geom_text(aes(label = p_val), data = sig, size = 10) +
    # scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-0.45, 0.45)) +
    scale_fill_gradient2(low = "blue", high = "red") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
