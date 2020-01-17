library(tidyverse)
library(mashr)


# Load the GWAS results ---------------------------------------------------
gwas <- list.files("data/gemma/output", "*\\.assoc\\.txt", full.names = TRUE) %>%
  map(read_tsv)
phenotypes <- list.files("data/gemma/output", "*\\.assoc\\.txt") %>%
  str_remove(., "\\.assoc\\.txt") %>%
  str_remove(., "norm_")

# Effect sizes
B <- map(gwas, function(df) pull(df, beta)) %>%
  bind_cols() %>% as.matrix()

# Standard errors of effect sizes
S <- map(gwas, function(df) pull(df, se)) %>%
  bind_cols() %>% as.matrix()

# Degrees of freedom
df <- 727 - 8 - 1

# Add phenotypes to the data matrices
names(gwas) <- colnames(B) <- colnames(S) <- phenotypes
rownames(B) <- rownames(S) <- gwas[[1]]$rs
write_rds(list(beta_hat = B, s_hat = S, df = df), "data/gemma/mash_data.rds")


# Initial mash analysis ---------------------------------------------------
# Identify the strongest effects and a random subset of effects to approximate
# the null distribution.
m_1by1 <- mash_1by1(mash_set_data(B, S, df = df))
write_rds(m_1by1, "data/gemma/mash_1by1.rds")
strong_subset <- get_significant_results(m_1by1, thresh = 0.1)
random_subset <- sample(1:nrow(B), 40000)

# Estimate the null correlation structure
data_temp <- mash_set_data(B[random_subset, ], S[random_subset, ], df = df)
Vhat <- estimate_null_correlation_simple(data_temp)
rm(data_temp); gc()

# Set up data structures
data_all <- mash_set_data(B, S, df = df, V = Vhat)
data_random <- mash_set_data(B[random_subset, ], S[random_subset, ], df = df, V = Vhat)
data_strong <- mash_set_data(B[strong_subset, ], S[strong_subset, ], df = df, V = Vhat)


# Derive data-driven covariance matrices ----------------------------------
# Column-centered Z-scores
Z <- apply(data_strong$Bhat/data_strong$Shat, 2, function(x) x - mean(x))

# Empirical covariance matrix
U1 <- (t(Z) %*% Z)/ncol(Z)

# Rank-P approximation of Z using SVD
P <- 3
svd_Z <- svd(Z, nu = P, nv = P)
U2 <- (svd_Z$v %*% diag(svd_Z$d[1:P])^2 %*% t(svd_Z$v))/ncol(Z)

# Rank-Q approximation of Z using SFA
Q <- 3
write_tsv(as_tibble(Z), "data/gemma/centered_cov.txt", col_names = FALSE)
system(paste0("./src/sfa/bin/sfa_linux -gen data/gemma/centered_cov.txt -g ", 
              nrow(Z), " -n ", ncol(Z), " -k ", Q, " -o data/gemma/gemma_norm"))

L <- read_tsv("data/gemma/gemma_norm_lambda.out", col_names = FALSE) %>% 
  as.matrix()
FA <- read_tsv("data/gemma/gemma_norm_F.out", col_names = FALSE) %>%
  as.matrix()

U3 <- (t(FA) %*% t(L) %*% L %*% FA)/ncol(Z)

# Extreme deconvolution of U1, U2, and U3
U_ed <- cov_ed(data_strong, list(U1, U2, U3))

# Rank-1 approximations of Z using SFA results
U_sfa <- lapply(1:Q, function(q) {
  t(FA[q, , drop = FALSE]) %*% t(L[, q, drop = FALSE]) %*%
    L[, q, drop = FALSE] %*% FA[q, , drop = FALSE]
})
names(U_sfa) <- paste0("SFA_", 1:Q)

# Canonical covariance matrices
U_c <- cov_canonical(data_random)


# Final mash analysis -----------------------------------------------------
m <- mash(data_random, Ulist = c(U_ed, U_sfa, U_c), outputlevel = 1)
write_rds(m, "data/gemma/mash_weights.rds")

# Compute posterior summaries for the strongest signals
m2 <- mash(data_strong, g = get_fitted_g(m), fixg = TRUE)
snp_names <- gwas[[1]]$rs[strong_subset]
for (i in 1:5) rownames(m2$result[[i]]) <- snp_names
write_rds(m2, "data/gemma/norm_snp_mash.rds")

# Compute posterior sumaries for all SNPs
m3 <- mash(data_all, g = get_fitted_g(m), fixg = TRUE)
for (i in 1:5) rownames(m3$result[[i]]) <- gwas[[1]]$rs
write_rds(m3, "data/gemma/norm_snp_mash_all.rds")
