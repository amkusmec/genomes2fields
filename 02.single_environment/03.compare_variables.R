library(tidyverse)
library(readxl)


# Compare variable inclusion ----------------------------------------------
vars1 <- read_xlsx("data/variable_inclusion.xlsx", na = "NA", skip = 1)
vars2 <- read_xlsx("data/variable_inclusion_gbs.xlsx", na = "NA", skip = 1)

# How many sites include each variable?
counts1 <- apply(vars1[, 3:8], 2, table)
counts2 <- apply(vars2[, 3:8], 2, table)

# Compare each site
similarity <- sapply(1:nrow(vars1), function(i) {
  all.equal(as.vector(t(vars1[i, ])), as.vector(t(vars2[i, ])))
})

table(similarity)
nrow(vars1) - sum(similarity == "TRUE")

# Compare each variable
var_sim <- sapply(3:8, function(i) sum(vars1[, i] == vars2[, i], na.rm = TRUE))
non_missing <- sapply(3:8, function(i) sum(!is.na(vars1[, i]) & !is.na(vars2[, i])))
non_missing - var_sim


# Compare hybrid BLUEs ----------------------------------------------------
res1 <- read_rds("data/phenotype/yield_stage_one.rds")
res2 <- read_rds("data/phenotype/yield_stage_one_gbs.rds")

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
ggsave("figures/single/blue_cor.pdf", width = 10, height = 5, units = "in", dpi = 300)

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
ggsave("figures/single/blue_rank.pdf", width = 10, height = 5, units = "in", dpi = 300)

abs_diff <- sapply(seq_along(res1), function(i) {
  temp <- inner_join(res1[[i]]$blue, res2[[i]]$blue, by = "PedigreeNew")
  median(abs(temp[[2]] - temp[[3]]))
})

tibble(Site = names(res1), MSD = abs_diff, Similarity = similarity) %>%
  arrange(MSD) %>%
  mutate(MSD = log10(MSD), 
         Site = factor(Site, levels = Site, ordered = TRUE), 
         Similarity = if_else(str_detect(Similarity, "1"), "1", 
                              if_else(str_detect(Similarity, "2"), "2", 
                                      if_else(str_detect(Similarity, "3"), "3", "0"))),
         Size = sapply(res2, function(x) nrow(x$blue))/
           sapply(res1, function(x) nrow(x$blue))) %>%
  ggplot(., aes(x = Site, y = MSD, colour = Similarity, shape = Similarity, size = Size)) + 
    theme_classic() + geom_point() + 
    labs(y = expression(paste(log[10], "(Median Absolute Difference)"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(type = "qual", palette = "Dark2")
ggsave("figures/single/blue_diff.pdf", width = 10, height = 5, units = "in", dpi = 300)


# Compare VCOV matrices ---------------------------------------------------
rs <- lapply(seq_along(res1), function(i) {
  hybrids <- rownames(res2[[i]]$vcov)
  temp <- res1[[i]]$vcov
  evolqg::RandomSkewers(temp[rownames(temp) %in% hybrids, colnames(temp) %in% hybrids], 
                        res2[[i]]$vcov, num.vectors = 1e4)
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
ggsave("figures/single/vcov_corr.pdf", width = 10, height = 5, units = "in", dpi = 300)

tibble(Site = names(res1), 
       MAD = log10(abs_diff), 
       Corr = rs$Corr) %>%
  ggplot(., aes(x = MAD, y = Corr)) + theme_classic() +
    geom_point(size = 2) +
    labs(x = expression(paste(log[10], "(Median Absolute Difference)")), 
         y = "Random Skewer Correlation")
ggsave("figures/single/mad_vs_corr.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Examine some particularly different sites -------------------------------
names(abs_diff) <- names(res1)
head(sort(abs_diff, decreasing = TRUE))

hybrids <- lapply(res2, function(x) rownames(x$vcov)) %>%
  unlist(use.names = FALSE) %>% unique()
yield <- read_rds("data/phenotype/yield_agron0.rds") %>%
  mutate(GBS = PedigreeNew %in% hybrids)

filter(yield, Site == "TXH2_2014") %>%
  ggplot(., aes(x = Column, y = Row, fill = Yield)) + theme_classic() +
    geom_tile(aes(size = GBS, colour = GBS)) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "yellow"))
ggsave("figures/single/TXH2_2014.pdf", width = 8, height = 5, units = "in", dpi = 300)
