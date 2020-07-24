### Not meant to be run directly.

### Implements a Bayesian random regression model on the selected environmental 
### variables.

library(tidyverse)
library(MCMCglmm)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-s", "--seed", type = "integer")
parser$add_argument("-c", "--chain", type = "integer")
args <- parser$parse_args()


# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015")

# g2 <- read_rds("data/weather/ga_final_model.rds")
sel <- read_lines("data/phenotype/selected_variables.txt")
data <- data[, c("BLUE", "PedigreeNew", "Site", sel)]

ped_site <- data %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# # Response
# y <- data$BLUE
# 
# # Predictors
# X <- model.matrix(as.formula(paste0("~ 1 + Site + ", 
#                                     paste(g2$g, collapse = " + "))), data = data)
# X[, -c(1:45)] <- scale(X[, -c(1:45)], center = TRUE, scale = TRUE)


# Variance-covariance matrix of the BLUEs
d <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
d <- d[names(d) %in% names(ped_site)]
d <- mapply(function(x, y) {
  idx <- which(rownames(x$vcov) %in% y$PedigreeNew)
  temp <- x$vcov[idx, idx]
  c(temp[1, 1], temp[1, 1] + diag(temp)[-1] + 2*temp[-1, 1])
}, x = d, y = ped_site)

W <- unlist(d, use.names = FALSE)
# data$SE <- sqrt(W)


# Random regression model -------------------------------------------------
prior <- list(R = list(V = diag(45)/2, n = 45), 
              G = list(G1 = list(V = diag(5)/2, n = 5)))

fixed_form <- paste0("BLUE ~ 1 + Site + ", 
                     paste(sel, collapse = " + ")) %>%
  as.formula()
rand_form <- paste0("~ us(1 + ", 
                    paste(sel, collapse = " + "), 
                    "):PedigreeNew") %>%
  as.formula()
rcov_form <- as.formula("~ idh(Site):units")

set.seed(args$seed)

start_time <- proc.time()[3]
m1 <- MCMCglmm(
  fixed = fixed_form, random = rand_form, rcov = rcov_form,
  prior = prior, pr = TRUE, family = "gaussian", mev = W,
  nitt = 2.5e6, burnin = 5e5, thin = 2e3, data = data, verbose = TRUE)
end_time <- proc.time()[3]

cat("Runtime:", (end_time - start_time)/60/60, "(h)")

write_rds(m1, paste0("data/phenotype/mcmc_chain", args$chain, ".rds"))
