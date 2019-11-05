library(tidyverse)
library(BGLR)
library(orthopolynom)
library(coda)
library(abind)

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


# Collect the output of the Gibbs sampler ---------------------------------
files <- list.files("data/bglr/rr_tmax", "*_TMAX[0-2]_b\\.bin", full.names = TRUE) %>% 
  split(., rep(1:4, each = 3))
mcmc_res_list <- lapply(files, function(f) { 
  temp <- lapply(f, function(f2) { 
      temp <- readBinMat(f2)
      if (nrow(temp) > 1000) temp <- temp[101:1100, ] 
      return(as.matrix(temp)) 
    })
  
  mcmc(do.call("cbind", temp), start = 1, end = 1000)
}) %>% 
  as.mcmc.list()


# Reformat the results into a tensor --------------------------------------
gibbs <- lapply(mcmc_res_list, as.matrix)
gibbs <- do.call("rbind", gibbs)
gibbs2 <- abind(list(gibbs[, 1:701], 
                     gibbs[, 702:1402], 
                     gibbs[, 1403:2103]), along = 3)
dimnames(gibbs2)[[3]] <- c("TMAX0", "TMAX1", "TMAX2")


# Calculate the posterior mean coefficient covariance matrix --------------
GC_list <- lapply(1:4000, function(i) cov(gibbs2[i, , ]))
GC_list <- do.call("abind", list(GC_list, along = 3))
GC <- apply(GC_list, c(1, 2), mean)


# Calculate the covariance function coefficients --------------------------
coefs <- polynomial.coefficients(leg_coef)
m <- max(sapply(coefs, length))
coefs <- lapply(coefs, function(f) c(f, rep(0, m - length(f))))
M <- do.call("rbind", coefs)
H <- t(M) %*% GC %*% M


# Calculate and plot the covariance function ------------------------------
tt <- seq(-1, 1, 0.01)
orig <- 0.5*(tt + 1)*(attr(data$TMAXS, "b") - attr(data$TMAXS, "a")) + attr(data$TMAXS, "a")
orig <- orig*sds[3] + means[3]
T <- matrix(c(rep(1, length(tt)), tt, tt^2), ncol = 3)
CF <- T %*% H %*% t(T)
CFc <- cov2cor(CF)
dCF <- as_tibble(cov2cor(CF), .name_repair = "unique") %>% 
  mutate(TMAXS1 = orig) %>% 
  gather(TMAXS2, CC, -TMAXS1) %>% 
  mutate(TMAXS2 = rep(orig, each = length(orig)))

ggplot(dCF, aes(x = TMAXS1, y = TMAXS2, fill = CC)) + theme_classic() + 
  geom_tile() + scale_fill_gradient2(midpoint = 0.5)

persp(z = CFc)
