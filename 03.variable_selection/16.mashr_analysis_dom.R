library(tidyverse)
library(mashr)

m2 <- read_rds("data/gemma/dominance/norm_dom_snp_mash.rds")


# Get the phenotypes ------------------------------------------------------
phenotypes <- list.files("data/gemma/dominance/output", "*\\.assoc\\.txt") %>%
  str_remove(., "norm_") %>%
  str_remove(., "_dom\\.assoc\\.txt")


# Mixture components for covariance structures ----------------------------
# Summarize mixture components
pi_mat <- matrix(m2$fitted_g$pi[-1], ncol = length(m2$fitted_g$Ulist), 
                 nrow = length(m2$fitted_g$grid), byrow = TRUE)
colnames(pi_mat) <- names(m2$fitted_g$Ulist)
colSums(pi_mat) %>% enframe(name = "Matrix") %>%
  bind_rows(., tibble(Matrix = "Null", value = m2$fitted_g$pi[1])) %>%
  ggplot(., aes(x = Matrix, y = value)) + theme_bw() +
    geom_col(colour = "black", fill = "red", alpha = 0.8) +
    labs(x = "", y = "Mixture proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.x = element_blank())
ggsave("figures/select/mash_mix_dom.pdf", width = 10, height = 5, units = "in", 
       dpi = 300)

# Summarize mixture components that account >0% of the components
colSums(pi_mat) %>% enframe(name = "Matrix") %>%
  filter(value > 0) %>%
  bind_rows(., tibble(Matrix = "Null", value = m2$fitted_g$pi[1])) %>%
  arrange(desc(value)) %>%
  mutate(Matrix = factor(Matrix, levels = Matrix, ordered = TRUE)) %>%
  ggplot(., aes(x = Matrix, y = value)) + theme_bw() +
    geom_col(colour = "black", fill = "red", alpha = 0.8) +
    labs(x = "", y = "Mixture proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.x = element_blank())
ggsave("figures/select/mash_mix0_dom.pdf", width = 8, height = 5, units = "in", 
       dpi = 300)


# Explore mixture components ----------------------------------------------
# Components that contribute >0% and are not single-phenotype components
comps <- c("ED_1", "ED_3", "ED_2")

for (cc in comps) {
  # Reform the covariance matrix
  x <- cov2cor(m2$fitted_g$Ulist[[cc]])
  rownames(x) <- colnames(x) <- colnames(m2$result$lfsr)
  x <- x[order(rownames(x)), order(colnames(x))]
  x[upper.tri(x, diag = FALSE)] <- NA
  
  # Correlation plot
  as_tibble(x, rownames = "Phenotype1") %>%
    gather(Phenotype2, R, -Phenotype1) %>%
    filter(!is.na(R)) %>%
    mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
           Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
    ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() +
      geom_tile() + labs(x = "", y = "", fill = "") +
      scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
  ggsave(paste0("figures/select/mash_", cc, "_dom.pdf"), width = 6, height = 4, 
         units = "in", dpi = 300)
  
  # Structure of the component
  vold <- svd(m2$fitted_g$Ulist[[cc]])$v
  u <- svd(m2$fitted_g$Ulist[[cc]])$u
  d <- svd(m2$fitted_g$Ulist[[cc]])$d
  v <- vold[order(colnames(m2$result$lfsr)), ]
  
  # Barplots of the eigenstructure
  for (j in 1:3) {
    tibble(Phenotype = rownames(x), 
           V = v[, j]/v[, j][which.max(abs(v[, j]))]) %>%
      ggplot(., aes(x = Phenotype, y = V)) + theme_classic() +
        geom_col(colour = "black") + labs(x = "", y = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste0(cc, ", Eigenvector ", j, ", (PVE = ", 
                       round(d[j]^2/sum(d^2), 2)*100, "%"))
    ggsave(paste0("figures/select/mash_", cc, "_dom_ev", j, ".pdf"), 
           width = 8, height = 5, units = "in", dpi = 300)
  }
}


# Number of significant phenotypes ----------------------------------------
thresh <- 0.001
n_sig <- get_n_significant_conditions(m2, thresh = thresh)
tibble(N_sig = n_sig) %>%
  filter(N_sig > 0) %>%
  ggplot(., aes(x = N_sig)) + theme_classic() +
    geom_histogram(binwidth = 1, colour = "black", fill = "aquamarine", 
                   alpha = 0.8) +
    labs(x = "Number of significant phenotypes", y = "Count")
ggsave("figures/select/sig_pheno_dom.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Overall sharing ---------------------------------------------------------
# Utility function from Urbut et al. (2019)
het.norm <- function(effectsize) {
  t(apply(effectsize, 1, function(x) {
    x/x[which.max(abs(x))]
  }))
}

m_data <- read_rds("data/gemma/dominance/mash_data.rds")
m_data[-3] <- lapply(m_data[-3], function(m) m[rownames(m2$result$lfsr), ])
pm_mash_beta <- m2$result$PosteriorMean*m_data$s_hat

# Sharing by sign
sig_mat <- m2$result$lfsr <= thresh
nsig <- rowSums(sig_mat)
sign_all <- mean(het.norm(pm_mash_beta[nsig > 0, ]) > 0) # 30.9%

# Sharing by magnitude
mag_all <- mean(het.norm(pm_mash_beta[nsig > 0, ]) > 0.5) # 14.3%


# Pairwise sharing by sign ------------------------------------------------
pm_mash_beta_mag <- pm_mash_beta[rowSums(m2$result$lfsr < thresh) > 0, ]
lfsr_mash <- m2$result$lfsr[rowSums(m2$result$lfsr < thresh) > 0, ]
shared_fold_size <- matrix(NA, nrow = ncol(lfsr_mash), ncol = ncol(lfsr_mash))
colnames(shared_fold_size) <- rownames(shared_fold_size) <- colnames(m2$result$lfsr)
for (i in 1:ncol(lfsr_mash)) {
  for (j in 1:ncol(lfsr_mash)) {
    sig_row <- which(lfsr_mash[, i] < thresh)
    sig_col <- which(lfsr_mash[, j] < thresh)
    a <- union(sig_row, sig_col)
    quotient <- pm_mash_beta_mag[a, i]/pm_mash_beta_mag[a, j]
    shared_fold_size[i, j] <- mean(quotient > 0)
  }
}
shared_fold_size <- shared_fold_size[order(rownames(shared_fold_size)), 
                                     order(colnames(shared_fold_size))]
shared_fold_size[upper.tri(shared_fold_size, diag = FALSE)] <- NA

as_tibble(shared_fold_size, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    scale_fill_distiller(type = "seq", palette = "Greens", 
                         limits = c(0, 1), direction = 1) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/select/mash_pairwise_sign_dom.pdf", width = 10, height = 8, 
       units = "in", dpi = 300)


# Pairwise sharing by magnitude -------------------------------------------
pm_mash_beta_mag <- pm_mash_beta[rowSums(m2$result$lfsr < thresh) > 0, ]
lfsr_mash <- m2$result$lfsr[rowSums(m2$result$lfsr < thresh) > 0, ]
shared_fold_size <- matrix(NA, nrow = ncol(lfsr_mash), ncol = ncol(lfsr_mash))
colnames(shared_fold_size) <- rownames(shared_fold_size) <- colnames(m2$result$lfsr)
for (i in 1:ncol(lfsr_mash)) {
  for (j in 1:ncol(lfsr_mash)) {
    sig_row <- which(lfsr_mash[, i] < thresh)
    sig_col <- which(lfsr_mash[, j] < thresh)
    a <- union(sig_row, sig_col)
    quotient <- pm_mash_beta_mag[a, i]/pm_mash_beta_mag[a, j]
    shared_fold_size[i, j] <- mean(quotient > 0.5 & quotient < 2)
  }
}
shared_fold_size <- shared_fold_size[order(rownames(shared_fold_size)), 
                                     order(colnames(shared_fold_size))]
shared_fold_size[upper.tri(shared_fold_size, diag = FALSE)] <- NA

as_tibble(shared_fold_size, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    scale_fill_distiller(type = "seq", palette = "Greens", 
                         limits = c(0, 1), direction = 1) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/select/mash_pairwise_mag_dom.pdf", width = 10, height = 8, 
       units = "in", dpi = 300)


# Phenotype-specific QTLs -------------------------------------------------
pm_mash_beta_norm <- het.norm(pm_mash_beta)
pm_mash_beta_norm <- pm_mash_beta[nsig > 0, ]
lfsr_mash <- m2$result$lfsr[nsig > 0, ]
a <- which(rowSums(pm_mash_beta_norm > 0.5) == 1)
lfsr_fold <- as.matrix(lfsr_mash[a, ])
pm <- as.matrix(pm_mash_beta_norm[a, ])
tspec <- NULL
for (i in 1:ncol(pm)) tspec[i] <- sum(pm[, i] > 0.5)
tspec <- as.matrix(tspec)
rownames(tspec) <- colnames(pm_mash_beta)


# Sharing summary ---------------------------------------------------------
# Another utility function
het.func <- function(normdat, threshold) {
  apply((normdat), 1, function(x) { sum(x > threshold) })
}

sign.func <- function(normeffectsize) {
  apply(normeffectsize, 1, function(x) sum(x > 0))
}

tibble(SNP = rownames(pm_mash_beta[nsig > 0, ]), 
       Magnitude = het.func(het.norm(pm_mash_beta[nsig > 0, ]), threshold = thresh), 
       Sign = sign.func(het.norm(pm_mash_beta[nsig > 0, ]))) %>%
  gather(Sharing, Count, Magnitude:Sign) %>%
  ggplot(., aes(x = Count, fill = Sharing)) + theme_classic() +
    geom_histogram(binwidth = 1, colour = "black", alpha = 0.8) +
    labs(x = "# Phenotypes", y = "# SNPs") +
    facet_wrap(~ Sharing) + guides(fill = "none") +
    scale_fill_manual(values = c("Magnitude" = "goldenrod", "Sign" = "limegreen"))
ggsave("figures/select/mash_pheno_sharing_dom.pdf", width = 8, height = 5, 
       units = "in", dpi = 300)


# Manhattan plot ----------------------------------------------------------
hlfsr <- apply(m2$result$lfsr, 1, min)
anno <- read_csv("data/gemma/snp_info.txt", col_names = FALSE) %>%
  rename(rs = X1, chr = X3, ps = X2) %>%
  full_join(., tibble(rs = names(hlfsr), pval = hlfsr), by = "rs") %>%
  mutate(pval = if_else(is.na(pval), 0.99999, pval))

source("src/manhattan_gemma.R")
p <- manhattan_gemma(anno, cutoff = -log10(thresh), ylab = expression(-log[10](min(lfsr))))
p
ggsave("figures/select/mash_manhattan_dom.pdf", width = 10, height = 5, 
       units = "in", dpi = 300)
