library(tidyverse)
library(mashr)


m2 <- read_rds("data/single_site_gwas/common_snp_mash.rds")


# Get the sites -----------------------------------------------------------
blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- names(blue)
sites <- sites[!str_detect(sites, "2017$")]
rm(blue); gc()


# Scrape information from the log files -----------------------------------
info <- seq_along(sites) %>%
  map_df(function(i) {
    txt <- read_lines(paste0("data/single_site_gwas/output/common_single", i, ".log.txt"))
    txt <- txt[str_detect(txt, "analyzed") | str_detect(txt, "pve")]
    txt <- str_replace(txt, "## [:print:]* = ", "")
    tibble(Site = sites[i], Samples = txt[1])
  }) %>%
  mutate(Samples = as.integer(Samples))


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
ggsave("figures/single/mash_mix.pdf", width = 10, height = 5, units = "in", dpi = 300)

# Summarize mixture components that account for > 0% of the components
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
ggsave("figures/single/mash_mix0.pdf", width = 8, height = 5, units = "in", dpi = 300)


# Explore major mixture components ----------------------------------------
# Components that contribute > 0% and are not single-site components
comps <- c("ED_1", "ED_2", "SFA_4")

for (cc in comps) {
  # Reformat the covariance matrix
  x <- cov2cor(m2$fitted_g$Ulist[[cc]])
  rownames(x) <- colnames(x) <- colnames(m2$result$lfsr)
  x <- x[order(rownames(x)), order(colnames(x))]
  x[upper.tri(x, diag = FALSE)] <- NA
  
  # Correlation plot
  as_tibble(x, rownames = "Site1") %>%
    gather(Site2, R, -Site1) %>%
    filter(!is.na(R)) %>%
    mutate(Site1 = factor(Site1, levels = sort(unique(Site1)), ordered = TRUE), 
           Site2 = factor(Site2, levels = rev(sort(unique(Site2))), ordered = TRUE)) %>%
    ggplot(., aes(x = Site1, y = Site2, fill = R)) + theme_classic() +
      geom_tile() + labs(x = "", y = "", fill = "") +
      scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
  ggsave(paste0("figures/single/mash_", cc, ".pdf"), width = 10, height = 8, 
         units = "in", dpi = 300)
  
  # Structure of the component
  vold <- svd(m2$fitted_g$Ulist[[cc]])$v
  u <- svd(m2$fitted_g$Ulist[[cc]])$u
  d <- svd(m2$fitted_g$Ulist[[cc]])$d
  v <- vold[order(colnames(m2$result$lfsr)), ]
  
  # Barplots of the eigenstructure
  for (j in 1:3) {
    tibble(Site = rownames(x), 
           V = v[, j]/v[, j][which.max(abs(v[, j]))]) %>%
      ggplot(., aes(x = Site, y = V)) + theme_classic() +
        geom_col(colour = "black") + labs(x = "", y = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste0(cc, ", Eigenvector ", j, ", (PVE = ", 
                       round(d[j]^2/sum(d^2), 2)*100, "%)"))
    ggsave(paste0("figures/single/mash_", cc, "_ev", j, ".pdf"), width = 10, 
           height = 5, units = "in", dpi = 300)
  }
}


# Number of significant conditions ----------------------------------------
n_sig <- get_n_significant_conditions(m2)
tibble(N_sig = n_sig) %>%
  filter(N_sig > 0) %>%
  ggplot(., aes(x = N_sig)) + theme_classic() +
    geom_histogram(binwidth = 1, colour = "black", fill = "aquamarine", 
                   alpha = 0.8) +
    labs(x = "Number of significant location-years", y = "Count")
ggsave("figures/single/sig_sites.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Effective sample sizes --------------------------------------------------
m_data <- read_rds("data/single_site_gwas/mash_data.rds")
m_data <- lapply(m_data, function(m) m[rownames(m2$result$lfsr), ])
ESS_factor <- apply(m_data$s_hat/m2$result$PosteriorSD, 2, median)
info <- tibble(Site = names(ESS_factor), ESS_factor = ESS_factor) %>%
  inner_join(., info, by = "Site") %>%
  mutate(ESS = Samples*ESS_factor, 
         Ratio = ESS/Samples) %>%
  select(Site, Samples, ESS_factor, ESS, Ratio)
info %>%
  select(-ESS_factor, -Ratio) %>%
  rename(`Sample Size` = Samples, `Effective Sample Size` = ESS) %>%
  gather(Type, Value, `Sample Size`, `Effective Sample Size`) %>%
  ggplot(., aes(x = Site, y = Value)) + theme_classic() +
    geom_col(aes(fill = Type), alpha = 0.8, colour = "black") +
    labs(x = "", y = "") + facet_wrap(~ Type) + guides(fill = "none") +
    coord_flip()
ggsave("figures/single/mash_ess.pdf", width = 10, height = 8, units = "in", dpi = 300)


# Overall sharing ---------------------------------------------------------
# Utility function from Urbut et al. (2019)
het.norm <- function(effectsize) {
  t(apply(effectsize,1,function(x){
    x/x[which.max(abs(x))]
  }))
}

pm_mash_beta <- m2$result$PosteriorMean*m_data$s_hat

# Sharing by sign
sig_mat <- m2$result$lfsr <= 0.05
nsig <- rowSums(sig_mat)
sign_all <- mean(het.norm(pm_mash_beta[nsig > 0, ]) > 0) # 77.8%

# Sharing by magnitude
mag_all <- mean(het.norm(pm_mash_beta[nsig > 0, ]) > 0.5) # 2.7%

### A follow-up to this would be to reassess sharing within subsets of the
### environments a la the brain/non-brain distinction in Urbut et al.


# Pairwise sharing by sign ------------------------------------------------
pm_mash_beta_mag <- pm_mash_beta[rowSums(m2$result$lfsr < 0.05) > 0, ]
lfsr_mash <- m2$result$lfsr[rowSums(m2$result$lfsr < 0.05) > 0, ]
shared_fold_size <- matrix(NA, nrow = ncol(lfsr_mash), ncol = ncol(lfsr_mash))
colnames(shared_fold_size) <- rownames(shared_fold_size) <- colnames(m2$result$lfsr)
for (i in 1:ncol(lfsr_mash)) {
  for (j in 1:ncol(lfsr_mash)) {
    sig_row <- which(lfsr_mash[, i] < 0.05)
    sig_col <- which(lfsr_mash[, j] < 0.05)
    a <- union(sig_row, sig_col)
    quotient <- pm_mash_beta_mag[a, i]/pm_mash_beta_mag[a, j]
    shared_fold_size[i, j] <- mean(quotient > 0)
  }
}
shared_fold_size <- shared_fold_size[order(rownames(shared_fold_size)), 
                                     order(colnames(shared_fold_size))]
shared_fold_size[upper.tri(shared_fold_size, diag = FALSE)] <- NA

as_tibble(shared_fold_size, rownames = "Site1") %>%
  gather(Site2, R, -Site1) %>%
  filter(!is.na(R)) %>%
  mutate(Site1 = factor(Site1, levels = sort(unique(Site1)), ordered = TRUE), 
         Site2 = factor(Site2, levels = rev(sort(unique(Site2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Site1, y = Site2, fill = R)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    scale_fill_distiller(type = "seq", palette = "Greens", 
                         limits = c(0, 1), direction = 1) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/single/mash_pairwise_sign.pdf", width = 10, height = 8, 
       units = "in", dpi = 300)


# Pairwise sharing by magnitude -------------------------------------------
pm_mash_beta_mag <- pm_mash_beta[rowSums(m2$result$lfsr < 0.05) > 0, ]
lfsr_mash <- m2$result$lfsr[rowSums(m2$result$lfsr < 0.05) > 0, ]
shared_fold_size <- matrix(NA, nrow = ncol(lfsr_mash), ncol = ncol(lfsr_mash))
colnames(shared_fold_size) <- rownames(shared_fold_size) <- colnames(m2$result$lfsr)
for (i in 1:ncol(lfsr_mash)) {
  for (j in 1:ncol(lfsr_mash)) {
    sig_row <- which(lfsr_mash[, i] < 0.05)
    sig_col <- which(lfsr_mash[, j] < 0.05)
    a <- union(sig_row, sig_col)
    quotient <- pm_mash_beta_mag[a, i]/pm_mash_beta_mag[a, j]
    shared_fold_size[i, j] <- mean(quotient > 0.5 & quotient < 2)
  }
}
shared_fold_size <- shared_fold_size[order(rownames(shared_fold_size)), 
                                     order(colnames(shared_fold_size))]
shared_fold_size[upper.tri(shared_fold_size, diag = FALSE)] <- NA

as_tibble(shared_fold_size, rownames = "Site1") %>%
  gather(Site2, R, -Site1) %>%
  filter(!is.na(R)) %>%
  mutate(Site1 = factor(Site1, levels = sort(unique(Site1)), ordered = TRUE), 
         Site2 = factor(Site2, levels = rev(sort(unique(Site2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Site1, y = Site2, fill = R)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    scale_fill_distiller(type = "seq", palette = "Greens", 
                         limits = c(0, 1), direction = 1) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/single/mash_pairwise_mag.pdf", width = 10, height = 8, 
       units = "in", dpi = 300)


# Site-specific QTLs ------------------------------------------------------
pm_mash_beta_norm <- het.norm(pm_mash_beta)
pm_mash_beta_norm <- pm_mash_beta_norm[nsig > 0, ]
lfsr_mash <- m2$result$lfsr[nsig > 0, ]
a <- which(rowSums(pm_mash_beta_norm > 0.5) == 1)
lfsr_fold <- as.matrix(lfsr_mash[a, ])
pm <- as.matrix(pm_mash_beta_norm[a, ])
tspec <- NULL
for (i in 1:ncol(pm)) tspec[i] <- sum(pm[, i] > 0.5)
tspec <- as.matrix(tspec)
rownames(tspec) <- colnames(pm_mash_beta)

as_tibble(tspec, rownames = "Site") %>%
  ggplot(., aes(x = Site, y = V1)) + theme_classic() +
    geom_col(colour = "black", fill = "steelblue", alpha = 0.8) +
    labs(x = "", y = "# Site-specific QTLs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/single/mash_specific_qtl.pdf", width = 10, height = 6, units = "in", dpi = 300)



# Sharing summary ---------------------------------------------------------
# Another utility function
het.func = function (normdat, threshold) {
  apply((normdat),1,function(x){sum(x > threshold)})
}

sign.func <- function (normeffectsize) {
  apply(normeffectsize,1,function(x)(sum(x>0)))
}

tibble(SNP = rownames(pm_mash_beta[nsig > 0, ]), 
       Magnitude = het.func(het.norm(pm_mash_beta[nsig > 0, ]), threshold = 0.5), 
       Sign = sign.func(het.norm(pm_mash_beta[nsig > 0, ]))) %>%
  gather(Sharing, Count, Magnitude:Sign) %>%
  ggplot(., aes(x = Count, fill = Sharing)) + theme_classic() +
    geom_histogram(binwidth = 1, colour = "black", alpha = 0.8) +
    labs(x = "# Sites", y = "# SNPs") +
    facet_wrap(~ Sharing) + guides(fill = "none") +
    scale_fill_manual(values = c("Magnitude" = "goldenrod", "Sign" = "limegreen"))
ggsave("figures/single/mash_site_sharing.pdf", width = 8, height = 5, 
       units = "in", dpi = 300)


# Manhattan plot ----------------------------------------------------------
# hlfsr <- apply(m2$result$lfsr, 1, function(x) length(x)/sum(1/x))
hlfsr <- apply(m2$result$lfsr[, -19], 1, min)
anno <- read_csv("data/single_site_gwas/snp_info.txt", col_names = FALSE) %>%
  rename(rs = X1, chr = X3, ps = X2) %>%
  full_join(., tibble(rs = names(hlfsr), pval = hlfsr), by = "rs") %>%
  mutate(pval = if_else(is.na(pval), 0.99999, pval))

source("src/manhattan_gemma.R")
p <- manhattan_gemma(anno, cutoff = -log10(0.05), ylab = expression(-log[10](ring(lfsr))))
p
ggsave("figures/single/mash_manhattan.pdf", width = 10, height = 5,
       units = "in", dpi = 300)
