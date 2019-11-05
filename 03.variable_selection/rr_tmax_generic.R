library(argparse)
library(tidyverse)
library(BGLR)
library(orthopolynom)
library(parallel)

parser <- ArgumentParser()
parser$add_argument("-d", "--degree", type = "integer")
args <- parser$parse_args()
degree <- args$degree

# Load the data -----------------------------------------------------------
# Phenotype and weather data
data <- read_rds("data/phenotype/yield_blue_final.rds") %>%
  separate(Site, c("Location", "Year"), sep = "_", remove = FALSE) %>%
  mutate(TMAX_X0.775_X0.875 = TMAX_X0.775_X0.875 - mean(TMAX_X0.775_X0.875))

ga <- c("NET_X0.025_X0.45", "NET_X0.65_X1.275", "TMAX_X0.775_X0.875", 
        "TMIN_X0.05_X1.5", "TMIN_X0.15_X0.65")


# Construct the model object ----------------------------------------------
data$TMAXS <- scaleX(data$TMAX_X0.775_X0.875, u = -1, v = 1)
leg_coef <- legendre.polynomials(n = degree, normalized = TRUE)
leg <- polynomial.values(leg_coef, data$TMAXS) %>%
  as.data.frame() %>% as_tibble()
names(leg) <- paste0("TMAX", 0:degree)

data <- bind_cols(data, leg)

if (degree < 2) {
  fix_fm <- "~ factor(Site)" %>% as.formula()
} else {
  fix_fm <- paste0("~ factor(Site) + ", 
                   paste(paste0("TMAX", 1:(degree - 1)), collapse = " + ")) %>%
    as.formula()
}

ETA <- list(fixed = list(fix_fm, data = data, model = "FIXED", saveEffects = TRUE))
for (i in 0:(degree - 1)) {
  ETA[[paste0("TMAX", i)]] <- list(paste0("~ 0 + factor(PedigreeNew):TMAX", i) %>%
                                    as.formula(), 
                                  data = data, model = "BRR", saveEffects = TRUE)
}

# Construct groups for heterogeneous error variances
groups <- as.integer(factor(data$Site))


# Compute the random regression model -------------------------------------
# Run 4 chains
cl <- makeCluster(4)
clusterEvalQ(cl, library(BGLR))
clusterExport(cl, list("ETA", "groups", "data", "degree"))

rxn <- parLapply(cl, 1:4, function(i) {
  BGLR(y = data$BLUE, ETA = ETA, burnIn = 1e4, nIter = 1.1e5, thin = 1e2, 
       saveAt = paste0("data/bglr/rr_tmax", degree, "/rr_tmax", 
                       degree, "_chain", i, "_"), groups = groups)
})

stopCluster(cl)

write_rds(rxn, paste0("data/bglr/rr_tmax", degree, "/rr_tmax", degree, "_models.rds"))
