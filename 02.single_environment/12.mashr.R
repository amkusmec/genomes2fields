library(tidyverse)
library(mashr)


# Load the GWAS results ---------------------------------------------------
gwas <- list.files("data/gemma/output", "^common_single[0-9]{1,2}\\.assoc\\.txt", 
                   full.names = TRUE) %>%
  map(read_tsv)
site_idx <- list.files("data/gemma/output", "^common_single[0-9]{1,2}\\.assoc\\.txt") %>%
  str_replace(., "common_single", "") %>%
  str_replace(., "\\.assoc\\.txt", "") %>%
  as.integer()

# Effect sizes
B <- map(gwas, function(df) pull(df, beta)) %>%
  bind_cols() %>% as.matrix()

# Standard errors of effect sizes
S <- map(gwas, function(df) pull(df, se)) %>%
  bind_cols() %>% as.matrix()

# Location-year-specific degrees of freedom
df <- list.files("data/gemma/output", "^common_single[0-9]{1,2}\\.log\\.txt", 
                 full.names = TRUE) %>%
  map_int(function(f) {
    temp <- read_lines(f)
    temp[str_detect(temp, "analyzed individuals")] %>%
      str_split(., " ") %>%
      sapply(., function(x) x[7]) %>%
      as.integer()
  })
df <- df - 9 - 1 # 9 principal components and the mean
df <- matrix(df, ncol = ncol(B), nrow = nrow(B), byrow = TRUE)

# Add location-year names to the data matrices
sites <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  pull(Site) %>%
  unique()
sites <- sites[site_idx]

names(gwas) <- colnames(B) <- colnames(S) <- colnames(df) <- sites
rownames(B) <- rownames(S) <- rownames(df) <- gwas[[1]]$rs
write_rds(list(beta_hat = B, s_hat = S, df = df), "data/gemma/mash_data.rds")


# Identify the strongest effects and a random subset of effects to approximate
# the null distribution
m_1by1 <- mash_1by1(mash_set_data(B, S, df = df))
strong_subset <- get_significant_results(m_1by1, thresh = 0.05)
random_subset <- sample(1:nrow(B), 40000)

# Estimate the null correlation structure
data_temp <- mash_set_data(B[random_subset, ], S[random_subset, ], 
                           df = df[random_subset, ])
Vhat <- estimate_null_correlation_simple(data_temp)
rm(data_temp); gc()

# Set-up data structures
data_random <- mash_set_data(B[random_subset, ], S[random_subset, ], 
                             df = df[random_subset, ], V = Vhat)
data_strong <- mash_set_data(B[strong_subset, ], S[strong_subset, ], 
                             df = df[strong_subset, ], V = Vhat)

# Data-driven covariance matrices
# Column-centered Z-scores
Z <- apply(data_strong$Bhat/data_strong$Shat, 2, function(x) x - mean(x))

# Empirical covariance matrix
U1 <- (t(Z) %*% Z)/ncol(Z)

# Rank-P approximation of Z using SVD
P <- 3
svd_Z <- svd(Z, nu = P, nv = P)
U2 <- (svd_Z$v %*% diag(svd_Z$d[1:P])^2 %*% t(svd_Z$v))/ncol(Z)

# Rank-Q approximation of Z using SFA
Q <- 5
write_tsv(as_tibble(Z), "data/gemma/centered_cov.txt", col_names = FALSE)
system(paste0("./src/sfa/bin/sfa_linux -gen data/gemma/centered_cov.txt -g ", 
              nrow(Z), " -n ", ncol(Z), " -k ", Q, " -o data/gemma/gemma_common"))

L <- read_tsv("data/gemma/gemma_common_lambda.out", col_names = FALSE) %>%
  as.matrix()
FA <- read_tsv("data/gemma/gemma_common_F.out", col_names = FALSE) %>%
  as.matrix()

U3 <- (t(FA) %*% t(L) %*% L %*% FA)/ncol(Z)

# Extreme deconvolution of U1, U2, and U3
U_ed <- cov_ed(data_strong, list(U1, U2, U3))

# Rank-1 approximations of Z using SFA results
U_sfa <- lapply(1:Q, function(q) t(FA[q, , drop = FALSE]) %*% t(L[, q, drop = FALSE]) %*% L[, q, drop = FALSE] %*% FA[q, , drop = FALSE])
names(U_sfa) <- paste0("SFA_", 1:Q)

# Canonical covariance matrices
U_c <- cov_canonical(data_random)

# Estimate mixture proportions
m <- mash(data_random, Ulist = c(U_ed, U_sfa, U_c), outputlevel = 1)

# Compute posterior summaries
m2 <- mash(data_strong, g = get_fitted_g(m), fixg = TRUE)
snp_names <- gwas[[1]]$rs[strong_subset]
for (i in 1:5) rownames(m2$result[[i]]) <- snp_names
write_rds(m2, "data/gemma/common_snp_mash.rds")
