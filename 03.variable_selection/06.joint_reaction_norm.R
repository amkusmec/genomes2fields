library(tidyverse)
library(BGLR)
library(parallel)
library(coda)


# Load the data -----------------------------------------------------------
# Phenotype and weather data
data <- read_rds("data/phenotype/yield_blue_final.rds") %>%
  separate(Site, c("Location", "Year"), sep = "_", remove = FALSE)

ga <- c("NET_X0.025_X0.45", "NET_X0.65_X1.275", "TMAX_X0.775_X0.875", 
        "TMIN_X0.05_X1.5", "TMIN_X0.15_X0.65")

# # Get the selected environmental variables
# ga <- read_rds("data/weather/ga_select_R.rds")[["TXH1_2016"]]$ga$g[-1]
# 
# tibble(V = ga) %>%
#   separate(V, c("Variable", "Start", "End"), sep = "_", remove = TRUE) %>%
#   mutate_at(c("Start", "End"), str_replace, pattern = "X", replacement = "") %>%
#   mutate_at(c("Start", "End"), as.numeric) %>%
#   mutate(Ypos = factor(Variable) %>% as.integer(), 
#          Ypos = Ypos + c(0.125, -0.125, 0, 0, 0)) %>%
#   ggplot(.) + theme_classic() +
#     geom_segment(aes(x = Start, xend = End, y = Ypos, yend = Ypos, 
#                      colour = Variable), size = 2) +
#     geom_vline(xintercept = 1, linetype = 2) +
#     labs(x = "% CHU to anthesis", y = "") + guides(colour = "none") +
#     scale_colour_manual(values = c("TMIN" = "blue", "TMAX" = "red", "PPT" = 
#                                      "skyblue", "NET" = "orange")) +
#     scale_y_continuous(breaks = 1:4, labels = c("NET", "PPT", "TMAX", "TMIN")) +
#     scale_x_continuous(limits = c(0, 1.5), labels = scales::percent)
# ggsave("figures/select/final_variables.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Center and scale the predictors -----------------------------------------
means <- apply(data[, 7:11], 2, mean)
sds <- apply(data[, 7:11], 2, sd)

for (i in 7:11) {
  data[[i]] <- (data[[i]] - means[i - 6])/sds[i - 6]
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
       saveAt = paste0("data/bglr/rxn_chain", i, "_"), groups = groups)
})

stopCluster(cl)

write_rds(rxn, "data/phenotype/joint_rxn_norm_model.rds")


# Check chain convergence -------------------------------------------------
site_levels <- sort(unique(data$Site))
files <- list.files("data/bglr", "*\\.", full.names = TRUE) %>%
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
        colnames(temp) <- str_replace(f2, "data/bglr/rxn_chain[1-4]_ETA_", "") %>%
          str_replace(., "\\.dat", "")
      } else if (str_detect(f2, "mu")) {
        colnames(temp) <- "mu"
      } else if (str_detect(f2, "varE")) {
        colnames(temp) <- site_levels[as.integer(names(rxn[[1]]$varE))]
      } else {
        comp <- str_replace(f2, "data/bglr/rxn_chain[1-4]_ETA_", "") %>%
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
length(reduce(idx, intersect)) # 3,377/4,309 (78%) parameters converged in all chains 
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

ggplot(ped, aes(x = Value)) + theme_classic() +
  geom_density() + facet_wrap(~ Variable, scales = "free") + 
  labs(x = "", y = "Density")
ggsave("figures/select/rxn_norm_density.pdf", width = 8, height = 5, units = "in", dpi = 300)

rxn_cor <- ped %>%
  spread(Variable, Value) %>%
  select(-PedigreeNew) %>%
  cor()
rxn_cor[upper.tri(rxn_cor, diag = FALSE)] <- NA
as_tibble(rxn_cor, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/select/rxn_norm_cor.pdf", width = 7, height = 5, units = "in", dpi = 300)

write_rds(ped, "data/phenotype/rxn_norm_parameters.rds")


# Estimate the effect of indices on environmental quality -----------------
e_main <- params %>%
  filter(str_detect(name, "Site")) %>%
  mutate(name = str_remove(name, "factor\\(Site\\)")) %>%
  rename(Site = name, Value = value) %>%
  bind_rows(., tibble(Site = "ARH1_2016", Value = 0)) %>%
  arrange(Site)
e_sum <- data %>%
  group_by(Site) %>%
  summarise_at(vars(NET_X0.025_X0.45:TMIN_X0.15_X0.65), median) %>%
  ungroup()
e_main <- inner_join(e_main, e_sum, by = "Site")

model <- lm(Value ~ . - 1 - Site, data = e_main)
enframe(coef(model)) %>%
  rename(Variable = name, Coefficient = value) %>%
  write_csv(., "data/phenotype/env_coef.csv")
