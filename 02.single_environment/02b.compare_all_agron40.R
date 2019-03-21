library(tidyverse)
library(readxl)


# Compare variable inclusion ----------------------------------------------
vars1 <- read_xlsx("data/variable_inclusion_all_agron0.xlsx", na = "NA", skip = 1)
vars2 <- read_xlsx("data/variable_inclusion_all_agron40.xlsx", na = "NA", skip = 1)

# How many sites include each variable?
counts1 <- apply(vars1[, 3:8], 2, table)
counts2 <- apply(vars2[, 3:8], 2, table)

# How many sites with variable changes?
na1 <- apply(vars1[, 3:8], 1, function(x) sum(is.na(x)))
na2 <- apply(vars2[, 3:8], 1, function(x) sum(is.na(x)))

sum(na1 != na2)

# How many variables were imputed per site?
table(na1 - na2)

# How many times was each variable imputed?
imputed <- sapply(3:5, function(i) {
  sum(is.na(vars1[[i]]) != is.na(vars2[[i]]))
})
imputed

# Similarities accounting for changes in missing variables
similarity <- sapply(which(na1 != na2), function(i) {
  all.equal(as.vector(t(vars1[i, 6:8])), as.vector(t(vars2[i, 6:8])))
})
table(similarity)

idx <- which(na1 != na2)[similarity != "TRUE"]
vars1[idx, ]; vars2[idx, ]

# For added variables, how many were included?
added <- sapply(3:5, function(i) {
  sum(vars2[[i]][is.na(vars1[[i]]) != is.na(vars2[[i]])] == "+")
})
added


# Compare hybrid BLUEs ----------------------------------------------------
res1 <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
res2 <- read_rds("data/phenotype/yield_stage_one_all_agron40.rds")

blue_cor <- sapply(seq_along(res1), function(i) {
  temp <- inner_join(res1[[i]]$blue, res2[[i]]$blue, by = "PedigreeNew")
  cor(temp[[2]], temp[[3]])
})

tibble(Site = names(res1), Pearson = blue_cor) %>%
  arrange(Pearson) %>%
  mutate(Site = factor(Site, levels = Site, ordered = TRUE)) %>%
  ggplot(., aes(x = Site, y = Pearson)) + theme_classic() +
    geom_point(size = 3) + labs(y = "r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/single/blue_cor_all_40.pdf", width = 10, height = 5, 
       units = "in", dpi = 300)

blue_rank <- sapply(seq_along(res1), function(i) {
  temp <- inner_join(res1[[i]]$blue, res2[[i]]$blue, by = "PedigreeNew")
  cor(temp[[2]], temp[[3]], method = "kendall")
})

tibble(Site = names(res1), Kendall = blue_rank) %>%
  arrange(Kendall) %>%
  mutate(Site = factor(Site, levels = Site, ordered = TRUE)) %>%
  ggplot(., aes(x = Site, y = Kendall)) + theme_classic() +
    geom_point(size = 3) + labs(y = expression(tau)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/single/blue_rank_all_40.pdf", width = 10, height = 5, 
       units = "in", dpi = 300)

abs_diff <- sapply(seq_along(res1), function(i) {
  temp <- inner_join(res1[[i]]$blue, res2[[i]]$blue, by = "PedigreeNew")
  median(abs(temp[[2]] - temp[[3]]))
})

tibble(Site = names(res1), MSD = abs_diff, Imputed = factor(na1 - na2)) %>%
  arrange(MSD) %>%
  mutate(MSD = log10(MSD + 1e-3), 
         Site = factor(Site, levels = Site, ordered = TRUE)) %>%
  ggplot(., aes(x = Site, y = MSD, colour = Imputed, size = Imputed)) +
    theme_classic() + geom_point() + geom_hline(yintercept = 0, linetype = 2) +
    labs(y = expression(paste(log[10], "(Median Absolute Difference)"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(type = "qual", palette = "Dark2")
ggsave("figures/single/blue_diff_all_40.pdf", width = 10, height = 5, 
       units = "in", dpi = 300)


# Compare VCOV matrices ---------------------------------------------------
rs <- lapply(seq_along(res1), function(i) {
  evolqg::RandomSkewers(res1[[i]]$vcov, res2[[i]]$vcov, num.vectors = 1e4)
})
rs <- tibble(Site = names(res1), 
             Corr = sapply(rs, function(x) x[1]), 
             P = sapply(rs, function(x) x[2]), 
             SD = sapply(rs, function(x) x[3]))

arrange(rs, Corr) %>%
  mutate(Site = factor(Site, levels = Site, ordered = TRUE)) %>%
  ggplot(., aes(x = Site, y = Corr)) + theme_classic() +
    geom_pointrange(aes(ymin = Corr - SD, ymax = Corr + SD)) +
    labs(y = "Random Skewer Correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/single/vcov_corr_all_40.pdf", width = 10, height = 5, 
       units = "in", dpi = 300)

tibble(Site = names(res1), 
       MAD = log10(abs_diff), 
       Corr = rs$Corr) %>%
  ggplot(., aes(x = MAD, y = Corr)) + theme_classic() +
    geom_point(size = 2) +
    labs(x = expression(paste(log[10], "(Median Absolute Difference)")), 
         y = "Random Skewer Correlation")
ggsave("figures/single/mad_vs_corr_all_40.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)


# Some very different sites -----------------------------------------------
tibble(No = res1$NYH2_2014$blue$BLUE, 
       Yes = res2$NYH2_2014$blue$BLUE) %>%
  ggplot(., aes(x = No, y = Yes)) + theme_classic() +
    geom_point(alpha = 0.8) + labs(x = "No Imputation", y = "Imputation") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
    ggtitle("NYH2_2014")

tibble(No = res1$MOH2_2014$blue$BLUE, 
       Yes = res2$MOH2_2014$blue$BLUE) %>%
  ggplot(., aes(x = No, y = Yes)) + theme_classic() +
    geom_point(alpha = 0.8) + labs(x = "No Imputation", y = "Imputation") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
    ggtitle("MOH2_2014")

# These are all instances where variables were imputed and then included in the 
# final model.
