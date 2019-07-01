library(tidyverse)
library(parallel)
library(purrrlyr)
library(QGenTools)


# Fit a monotonic decreasing spline to the estimated LD -------------------
# Code from:
# https://gist.github.com/jnpaulson/c47f9bd3246f1121ad3a911e5f707355
est_distance <- function(ld, threshold = 0.3, k = 30) {
  require(mgcv)
  
  init_gam <- gam(R2 ~ s(Bin, k = k, bs = "cr"), data = ld)
  sm <- smoothCon(s(Bin, k = k, bs = "cr"), ld, knots = NULL)[[1]]
  mc <- mono.con(sm$xp, up = FALSE)
  M <- list(X = sm$X, y = ld$R2, C = matrix(0, 0, 0), Ain = mc$A, 
            bin = mc$b, sp = init_gam$sp, p = -sm$xp, S = sm$S, 
            w = ld$R2*0 + 1, off = 0)
  p <- pcls(M)
  which.min(abs(drop(Predict.matrix(sm, data.frame(Bin = ld$Bin)) %*% p) - threshold))
}


# Load the data -----------------------------------------------------------
snps <- read_rds("data/gbs/add_snps.rds")
struct <- read_rds("data/gbs/pca_covariates.rds")


# Some parameters ---------------------------------------------------------
window_size <- 400 # In SNPs
step_size   <-  25 # In SNPs
bin_size    <- 2e4 # In kbp
AB <- realized_ab(snps$GD)


# Calculate corrected LD --------------------------------------------------
cl <- makeCluster(10)
clusterEvalQ(cl, { library(tidyverse); library(purrrlyr); library(QGenTools) })
clusterExport(cl, list("step_size", "window_size", "bin_size", "snps", "AB", "struct"))
ld <- parLapply(cl, split(snps$GM, snps$GM$Chromosome), function(df) {
    tibble(Start = seq(1, (nrow(df) %/% step_size)*step_size - window_size + 1, step_size), 
           End = seq(window_size, (nrow(df) %/% step_size)*step_size, step_size)) %>%
      by_row(function(r) {
        r2vs(snps$GD[, colnames(snps$GD) %in% df$SNP[r$Start[1]:r$End[1]]], AB, struct) %>%
          separate(Locus1, c("Chr1", "Pos1"), sep = "_") %>%
          separate(Locus2, c("Chr2", "Pos2"), sep = "_") %>%
          mutate_at(c("Pos1", "Pos2"), as.integer) %>%
          select(-Chr1, -Chr2)
      }) %>%
      unnest(.out) %>%
      mutate(Distance = (Pos2 - Pos1)/bin_size, 
             Bin = floor(Distance)) %>%
      group_by(Bin) %>%
      summarise(R2 = quantile(R2, 0.95)) %>%
      ungroup() %>%
      mutate(Chromosome = df$Chromosome[1]) %>%
      select(Chromosome, everything())
  })
ld <- bind_rows(ld)
write_rds(ld, "data/gbs/ld_decay.rds")

# Plot the output for a visual sanity check
ld %>%
  mutate(Bin = Bin*bin_size/1e6) %>%
  ggplot(., aes(x = Bin, y = R2)) + theme_classic() +
    geom_point(colour = "grey", size = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0.2, linetype = 2, colour = "red") +
    labs(x = "Position (Mbp)", y = expression(r[VS]^2)) +
    facet_wrap(~ Chromosome)


# Calculate uncorrected LD ------------------------------------------------
ld_unc <- parLapply(cl, split(snps$GM, snps$GM$Chromosome), function(df) {
    tibble(Start = seq(1, (nrow(df) %/% step_size)*step_size - window_size + 1, step_size), 
           End = seq(window_size, (nrow(df) %/% step_size)*step_size, step_size)) %>%
      by_row(function(r) {
        r2(snps$GD[, colnames(snps$GD) %in% df$SNP[r$Start[1]:r$End[1]]]) %>%
          separate(Locus1, c("Chr1", "Pos1"), sep = "_") %>%
          separate(Locus2, c("Chr2", "Pos2"), sep = "_") %>%
          mutate_at(c("Pos1", "Pos2"), as.integer) %>%
          select(-Chr1, -Chr2)
      }) %>%
      unnest(.out) %>%
      mutate(Distance = (Pos2 - Pos1)/bin_size, 
             Bin = floor(Distance)) %>%
      group_by(Bin) %>%
      summarise(R2 = quantile(R2, 0.95)) %>%
      ungroup() %>%
      mutate(Chromosome = df$Chromosome[1]) %>%
      select(Chromosome, everything())
  })
ld_unc <- bind_rows(ld_unc)
write_rds(ld_unc, "data/gbs/ld_unc_decay.rds")
stopCluster(cl)

# Plot the output for a visual sanity check
ld_unc %>%
  mutate(Bin = Bin*bin_size/1e6) %>%
  ggplot(., aes(x = Bin, y = R2)) + theme_classic() +
    geom_point(colour = "grey", size = 0.5, alpha = 0.8) +
    labs(x = "Position (Mbp)", y = expression(r^2)) +
    facet_wrap(~ Chromosome)

# Plot the differences to assess the magnitude of the effect of structure
# correction
ld %>%
  mutate(R2unc = ld_unc$R2, 
         Diff = R2unc - R2, 
         Bin = Bin*bin_size/1e6) %>%
  ggplot(., aes(x = Bin, y = Diff)) + theme_classic() +
    geom_point(colour = "grey", size = 0.5, alpha = 0.8) + 
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    labs(x = "Position (Mbp)", y = expression(r^2-r[VS]^2)) +
    facet_wrap(~ Chromosome)

# Plot the pairwise estimates
ld %>%
  mutate(R2unc = ld_unc$R2) %>%
  ggplot(., aes(x = R2, y = R2unc)) + theme_classic() +
    geom_point(colour = "grey", size = 0.5, alpha = 0.8) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
    labs(x = expression(r[VS]^2), y = expression(r^2)) +
    facet_wrap(~ Chromosome)


# Estimate the chromosome-wide decay rates --------------------------------
decay <- ld %>%
  split(., .$Chromosome) %>%
  map_dbl(est_distance, threshold = 0.3)
decay <- decay*bin_size/1e6

decay_unc <- ld_unc %>%
  split(., .$Chromosome) %>%
  map_dbl(est_distance, threshold = 0.3)
decay_unc <- decay_unc*bin_size/1e6


# An example smoothing plot -----------------------------------------------
df <- filter(ld, Chromosome == 5)
init_gam <- gam(R2 ~ s(Bin, k = 30, bs = "cr"), data = df)
sm <- smoothCon(s(Bin, k = 30, bs = "cr"), df, knots = NULL)[[1]]
mc <- mono.con(sm$xp, up = FALSE)
M <- list(X = sm$X, y = df$R2, C = matrix(0, 0, 0), Ain = mc$A, 
          bin = mc$b, sp = init_gam$sp, p = -sm$xp, S = sm$S, 
          w = df$R2*0 + 1, off = 0)
p <- pcls(M)
pred <- drop(Predict.matrix(sm, data.frame(Bin = df$Bin)) %*% p)
cutoff <- which.min(abs(pred - 0.3))

pdf("figures/munge/decay_example.pdf", width = 6, height = 4)
visreg::visreg(init_gam, "Bin", main = "Chromosome 5")
lines(df$Bin, pred, lty = 2, lwd = 3, col = "yellow")
abline(h = 0.3, lty = 3, lwd = 2, col = "red")
abline(v = cutoff, lty = 3, lwd = 2, col = "red")
dev.off()
