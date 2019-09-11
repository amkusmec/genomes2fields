library(tidyverse)
library(BGLR)
library(orthopolynom)
library(parallel)
library(coda)


# Load the data -----------------------------------------------------------
# Phenotype and weather data
data <- read_rds("data/phenotype/yield_blue_final.rds") %>%
  separate(Site, c("Location", "Year"), sep = "_", remove = FALSE)

ga <- c("NET_X0.025_X0.45", "NET_X0.65_X1.275", "TMAX_X0.775_X0.875", 
        "TMIN_X0.05_X1.5", "TMIN_X0.15_X0.65")


# Center and scale the predictors -----------------------------------------
means <- apply(data[, 7:11], 2, mean)
sds <- apply(data[, 7:11], 2, sd)

for (i in 7:11) data[[i]] <- (data[[i]] - means[i - 6])/sds[i - 6]


# Construct the model object ----------------------------------------------
data$TMAXS <- scaleX(data$TMAX_X0.775_X0.875, u = -1, v = 1)
leg_coef <- legendre.polynomials(n = 2, normalized = TRUE)
leg <- polynomial.values(leg_coef, data$TMAXS) %>%
  as.data.frame() %>% as_tibble()
names(leg) <- c("TMAX0", "TMAX1", "TMAX2")

data <- bind_cols(data, leg)

ETA <- list(fixed = list(~ factor(Site) + TMAXS, data = data, model = "FIXED", 
                         saveEffects = TRUE), 
            TMAX0 = list(~ 0 + factor(PedigreeNew):TMAX0, data = data, 
                         model = "BRR", saveEffects = TRUE), 
            TMAX1 = list(~ 0 + factor(PedigreeNew):TMAX1, data = data, 
                         model = "BRR", saveEffects = TRUE), 
            TMAX2 = list(~ 0 + factor(PedigreeNew):TMAX2, data = data, 
                         model = "BRR", saveEffects = TRUE))

# Construct groups for heterogeneous error variances
groups <- as.integer(factor(data$Site))


# Compute the random regression model -------------------------------------
# Run 4 chains
cl <- makeCluster(4)
clusterEvalQ(cl, library(BGLR))
clusterExport(cl, list("ETA", "groups", "data"))

rxn <- parLapply(cl, 1:4, function(i) {
  BGLR(y = data$BLUE, ETA = ETA, burnIn = 1e4, nIter = 1.1e5, thin = 1e2, 
       saveAt = paste0("data/bglr/rr_tmax/rr_tmax_chain", i, "_"), groups = groups)
})

stopCluster(cl)


# Check chain convergence -------------------------------------------------
site_levels <- sort(unique(data$Site))
files <- list.files("data/bglr/rr_tmax", "*", full.names = TRUE) %>%
  split(., rep(1:4, each = 9))
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
        colnames(temp) <- str_remove(f2, "data/bglr/rr_tmax/rr_tmax_chain[1-4]_ETA_") %>%
          str_remove(., "\\.dat")
      } else if (str_detect(f2, "mu")) {
        colnames(temp) <- "mu"
      } else if (str_detect(f2, "varE")) {
        colnames(temp) <- site_levels[as.integer(names(rxn[[1]]$varE))]
      } else {
        comp <- str_remove(f2, "data/bglr/rr_tmax/rr_tmax_chain[1-4]_ETA_") %>%
          str_remove(., "_b\\.bin")
        colnames(temp) <- rxn[[1]]$ETA[[comp]]$colNames %>%
          str_remove(., "factor\\(PedigreeNew\\)")
      }
    }
    
    if (nrow(temp) > 1000) temp <- temp[101:1100, ]
    
    return(as.matrix(temp))
  })
  
  mcmc(do.call("cbind", temp), start = 1, end = 1000)
}) %>%
  as.mcmc.list()

# 2,199 parameters total
gdiag <- geweke.diag(mcmc_res_list)
idx <- lapply(gdiag, function(g) which(abs(g$z) < 1.96))
length(reduce(idx, intersect)) # 1,741 (79%) parameters converged in all chains
length(reduce(idx, union)) # 100% parameters converged in at least one chain
unlist(idx) %>% enframe() %>% count(value) %>% count(n)
# Nothing converges in only one chain


# Collect the Gibbs samples -----------------------------------------------
gibbs <- lapply(mcmc_res_list, as.matrix)
gibbs <- do.call("rbind", gibbs)
params <- apply(gibbs, 2, mean) %>%
  enframe()

mus <- params$value[params$name == "TMAXS"]
ped <- filter(params, str_detect(name, "/")) %>%
  separate("name", c("PedigreeNew", "Variable"), sep = ":", remove = TRUE) %>%
  mutate(Variable = if_else(is.na(Variable), "Hybrid", Variable)) %>%
  rename(Value = value)

# Compute MSE for each hybrid
mu <- mean(gibbs[, "mu"])
model_mat <- cbind(model.matrix(~ 0 + factor(Site) + TMAXS, data = data)[, -1], 
                   model.matrix(~ 0 + factor(PedigreeNew):TMAX0, data = data), 
                   model.matrix(~ 0 + factor(PedigreeNew):TMAX1, data = data), 
                   model.matrix(~ 0 + factor(PedigreeNew):TMAX2, data = data))
colnames(model_mat) <- if_else(str_detect(colnames(model_mat), "factor"), 
                               str_remove(colnames(model_mat), "factor\\(PedigreeNew\\)"), 
                               colnames(model_mat))

beta_hat <- params %>%
  slice(1:2151) %>%
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


plot_hyb <- function(hyb) {
  # Predicted reaction norm
  rr <- seq(-1, 1, length.out = 500)
  y_hat <- mu + mus*rr + 
    params$value[params$name == paste(hyb, "TMAX0", sep = ":")]*0.7071 + 
    params$value[params$name == paste(hyb, "TMAX1", sep = ":")]*(1.2247*rr) + 
    params$value[params$name == paste(hyb, "TMAX2", sep = ":")]*(2.3717*rr^2 - 0.7906)
  
  # Confidence interval
  d <- matrix(NA, nrow = 4000, ncol = 500)
  for (i in seq_along(rr)) {
    d[, i] <- mu + mus*rr[i] + gibbs[, paste(hyb, "TMAX0", sep = ":")]*0.7071 +
      gibbs[, paste(hyb, "TMAX1", sep = ":")]*(1.2247*rr[i]) +
      gibbs[, paste(hyb, "TMAX2", sep = ":")]*(2.3717*rr[i]^2 - 0.7906)
  }
  q <- apply(d, 2, quantile, probs = c(0.025, 0.975)) %>%
    t() %>% as_tibble() %>%
    mutate(TMAXS = rr)
  
  # Observed yields adjusted for environment effects
  idx <- which(data$PedigreeNew == hyb)
  idx2 <- match(paste0("factor(Site)", data$Site[idx]), params$name)
  y <- data$BLUE[idx] - params$value[idx2]
  y <- c(data$BLUE[idx[1]], data$BLUE[idx[-1]] - params$value[idx2[-1]])
  
  # Plot it
  ggplot() + theme_bw() + labs(x = "Scaled Max. Temp.", y = "BLUE") + 
    geom_ribbon(aes(x = TMAXS, ymin = `2.5%`, ymax = `97.5%`), 
                data = q, fill = "grey70", alpha = 0.5) +
    geom_point(aes(x = TMAXS, y = Y), 
               data = tibble(TMAXS = data$TMAXS[idx], Y = y)) +
    geom_line(aes(x = TMAXS, y = Y), 
              data = tibble(TMAXS = rr, Y = y_hat), linetype = 2)
}

hyb_count <- data %>% count(PedigreeNew) %>% arrange(desc(n))
plot_hyb("LH195/LH185")

evars <- filter(params, str_detect(name, "_"), !str_detect(name, "var"), 
                !str_detect(name, "Site"), !str_detect(name, "/"))
vars <- filter(params, str_detect(name, "var"))
evars$TMAX0 <- vars$value[1]
evars$TMAX1 <- vars$value[2]
evars$TMAX2 <- vars$value[3]

evars %>%
  rename(Site = name, Residual = value) %>%
  gather(Factor, Value, Residual:TMAX2) %>%
  group_by(Site) %>%
  mutate(Value = Value/sum(Value)) %>%
  ungroup() %>%
  mutate(Factor = factor(Factor, ordered = TRUE, 
                         levels = rev(c("TMAX0", "TMAX1", "TMAX2", "Residual")))) %>%
  ggplot(., aes(x = Site, y = Value)) + theme_classic() +
    geom_col(aes(fill = Factor), colour = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # scale_fill_brewer(type = "qual", palette = "Set2") +
    scale_y_continuous(labels = scales::percent)

summary(sum(vars$value)/(sum(vars$value) + evars$value))

ped %>% spread(Variable, Value) %>%
  select(-PedigreeNew) %>%
  cor()

zero <- ped %>%
  spread(Variable, Value) %>%
  mutate(Zero = -TMAX1/(2*TMAX2)) %>%
  filter(Zero >= -1, Zero <= 1) %>%
  mutate(Scaled = (2.71 - -2)*(Zero - -1)/(1 - -1) + -2, 
         TMAX = Scaled*sds[3] + means[3], 
         Dir = if_else(TMAX2 > 0, "Up", if_else(TMAX2 < 0, "Down", "???")))

ggplot(zero, aes(x = TMAX, colour = Dir)) + geom_density()
ggplot(zero, aes(x = TMAX0, colour = Dir)) + geom_density()
