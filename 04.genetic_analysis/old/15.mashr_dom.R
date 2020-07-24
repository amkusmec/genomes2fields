library(tidyverse)
library(mashr)


# Load the GWAS results ---------------------------------------------------
gwas <- list.files("data/gemma/dominance/output", "*\\.assoc\\.txt", full.names = TRUE) %>%
  map(read_tsv)
phenotypes <- list.files("data/gemma/dominance/output", "*\\.assoc\\.txt") %>%
  str_remove(., "\\.assoc\\.txt") %>%
  str_remove(., "norm_")

# Effect sizes
B <- map(gwas, function(df) pull(df, beta)) %>%
  bind_cols() %>% as.matrix()

# Standard errors of effect sizes
S <- map(gwas, function(df) pull(df, se)) %>%
  bind_cols() %>% as.matrix()

# Degrees of freedom
df <- 701 - 8 - 1

# Add phenotypes to the data matrices
names(gwas) <- colnames(B) <- colnames(S) <- phenotypes
rownames(B) <- rownames(S) <- gwas[[1]]$rs
write_rds(list(beta_hat = B, s_hat = S, df = df), "data/gemma/dominance/mash_data.rds")


# Initial mash analysis ---------------------------------------------------
# Identify the strongest effects and a random subset of effects to approximate
# the null distribution.
m_1by1 <- mash_1by1(mash_set_data(B, S, df = df))
strong_subset <- get_significant_results(m_1by1, thresh = 0.1)

# No SNPs meet the threshold for the one-by-one analysis
# Proceed using all data
data <- mash_set_data(B, S, df = df)


# Derive data-driven covariance matrices ----------------------------------
# Column-centered Z-scores
Z <- apply(data$Bhat/data$Shat, 2, function(x) x - mean(x))

# Empirical covariance matrix
U1 <- (t(Z) %*% Z)/ncol(Z)

# Rank-P approximation of Z using SVD
P <- 3
svd_Z <- svd(Z, nu = P, nv = P)
U2 <- (svd_Z$v %*% diag(svd_Z$d[1:P])^2 %*% t(svd_Z$v))/ncol(Z)

# Rank-Q approximation of Z using SFA and 20k random markers
idx <- sample(1:nrow(Z), 20000)
Q <- 5
write_tsv(as_tibble(Z[idx, ]), "data/gemma/dominance/centered_cov_dom.txt", col_names = FALSE)
system(paste0("./src/sfa/bin/sfa_linux -gen data/gemma/dominance/centered_cov_dom.txt -g ", 
              length(idx), " -n ", ncol(Z), " -k ", Q, " -o data/gemma/dominance/gemma_norm_dom"))

L <- read_tsv("data/gemma/dominance/gemma_norm_dom_lambda.out", col_names = FALSE) %>%
  as.matrix()
FA <- read_tsv("data/gemma/dominance/gemma_norm_dom_F.out", col_names = FALSE) %>%
  as.matrix()

U3 <- (t(FA) %*% t(L) %*% L %*% FA)/ncol(Z)

# Extreme deconvolution of U1, U2, and U3
U_ed <- cov_ed(data, list(U1, U2, U3))

# Rank-1 approximations of Z using SFA results
U_sfa <- lapply(1:Q, function(q) {
  t(FA[q, , drop = FALSE]) %*% t(L[, q, drop = FALSE]) %*%
    L[, q, drop = FALSE] %*% FA[q, , drop = FALSE]
})
names(U_sfa) <- paste0("SFA_", 1:Q)

# Canonical covariance matrices
U_c <- cov_canonical(data)


# Final mash analysis -----------------------------------------------------
m <- mash(data, Ulist = c(U_ed, U_sfa, U_c), outputlevel = 2)
for (i in 1:5) rownames(m$result[[i]]) <- gwas[[1]]$rs
write_rds(m, "data/gemma/dominance/norm_dom_snp_mash.rds")
