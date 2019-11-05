library(tidyverse)
library(BGLR)
library(orthopolynom)
library(coda)
library(abind)

# Not sure if this part is correct
# models <- paste0("data/bglr/rr_tmax", 1:5, "/rr_tmax", 1:5, "_models.rds") %>%
#   map(read_rds)
# names(models) <- paste0("TMAX", 1:5)
# ll <- sapply(models, function(l) {
#   sapply(l, function(x) x$fit$postMeanLogLik) %>% mean
# })
# df <- (46 + 1 + 45 + 701*(1:5 - 1) + (1:5 - 1))
# 
# pchisq(-2*(ll[1] - ll[3]), df[3] - df[1])
# pchisq(-2*(ll[3] - ll[5]), df[5] - df[3])
# pchisq(-2*(ll[2] - ll[4]), df[4] - df[2])


# Load the data -----------------------------------------------------------
# Phenotype and weather data
data <- read_rds("data/phenotype/yield_blue_final.rds") %>%
  separate(Site, c("Location", "Year"), sep = "_", remove = FALSE) %>%
  mutate(TMAX_X0.775_X0.875 = TMAX_X0.775_X0.875 - mean(TMAX_X0.775_X0.875))

ga <- c("NET_X0.025_X0.45", "NET_X0.65_X1.275", "TMAX_X0.775_X0.875", 
        "TMIN_X0.05_X1.5", "TMIN_X0.15_X0.65")

# Construct the model object ----------------------------------------------
data$TMAXS <- scaleX(data$TMAX_X0.775_X0.875, u = -1, v = 1)
leg_coef <- legendre.polynomials(n = 5, normalized = TRUE)
leg <- polynomial.values(leg_coef, data$TMAXS) %>%
  as.data.frame() %>% as_tibble()
names(leg) <- paste0("TMAX", 0:5)

data <- bind_cols(data, leg)


# Collect the output of the Gibbs sampler ---------------------------------
res <- lapply(1:5, function(i) {
  files <- list.files(paste0("data/bglr/rr_tmax", i), 
                      "*_TMAX[0-4]_b\\.bin", full.names = TRUE) %>%
    split(., rep(1:4, each = i))
  
  lapply(files, function(f) {
    temp <- lapply(f, function(f2) { 
        temp <- readBinMat(f2)
        if (nrow(temp) > 1000) temp <- temp[101:1100, ]
        return(as.matrix(temp))
      })
    
    mcmc(do.call("cbind", temp), start = 1, end = 1000)
  }) %>%
    as.mcmc.list()
})
names(res) <- paste0("TMAX", 0:4)


# Reformat the results into a tensor --------------------------------------
gibbs <- lapply(res, function(mc) {
  do.call("rbind", lapply(mc, as.matrix))
})

gibbs2 <- lapply(gibbs, function(g) {
  k <- ncol(g)/701
  g2 <- lapply(1:k, function(i) g[, ((i - 1)*701 + 1):(i*701)])
  abind(g2, along = 3)
})

gibbs2 <- lapply(gibbs2, function(g) {
  dimnames(g)[[3]] <- paste0("TMAX", 0:(dim(g)[3] - 1))
  g
})


# Calculate the posterior mean coefficient covariance matrices ------------
GC_list <- lapply(gibbs2, function(g) {
  if (dim(g)[3] > 1) {
    do.call("abind", list(lapply(1:4000, function(i) cov(g[i, , ])), along = 3))
  } else {
    apply(g, 1, var)
  }
})

GC <- lapply(GC_list, function(g) {
  if (is.null(dim(g))) {
    mean(g)
  } else {
    apply(g, c(1, 2), mean)
  }
})


# Calculate the covariance function coefficients --------------------------
coefs <- polynomial.coefficients(leg_coef)
m <- max(sapply(coefs, length))
coefs <- lapply(coefs, function(f) c(f, rep(0, m - length(f))))
M <- do.call("rbind", coefs)
H_list <- lapply(GC, function(gc) {
  if (is.null(dim(gc))) {
    gc*M[1, 1]^2
  } else {
    m <- ncol(gc)
    t(M[1:m, 1:m]) %*% gc %*% M[1:m, 1:m]
  }
})


# Calculate the covariance functions --------------------------------------
tt <- seq(-1, 1, 0.01)
TT <- matrix(c(rep(1, length(tt)), tt, tt^2, tt^3, tt^4), ncol = 5)
CF_list <- lapply(H_list, function(h) {
  if (is.null(dim(h))) {
    h*(TT[, 1, drop = FALSE] %*% t(TT[, 1, drop = FALSE]))
  } else {
    m <- ncol(h)
    TT[, 1:m] %*% h %*% t(TT[, 1:m])
  }
})

orig <- 0.5*(tt + 1)*(attr(data$TMAXS, "b") - attr(data$TMAXS, "a")) + attr(data$TMAXS, "a")
dCF <- map_df(seq_along(CF_list), function(i) {
  as_tibble(cov2cor(CF_list[[i]]), .name_repair = "unique") %>%
    mutate(TMAXS1 = orig) %>%
    gather(TMAXS2, CC, -TMAXS1) %>%
    mutate(TMAXS2 = rep(orig, each = length(orig)), 
           K = paste0("k = ", i))
})


ggplot(dCF, aes(x = TMAXS1, y = TMAXS2, fill = CC)) + theme_classic() +
  geom_tile() + scale_fill_gradient2(midpoint = mean(range(dCF$CC))) + 
  facet_wrap(~ K)

filter(dCF, K == "k = 5") %>%
  ggplot(., aes(x = TMAXS1, y = TMAXS2, fill = CC)) +
    theme_classic() + geom_tile() + 
    scale_fill_gradient2(midpoint = mean(range(dCF$CC[dCF$K == "k = 5"])))
