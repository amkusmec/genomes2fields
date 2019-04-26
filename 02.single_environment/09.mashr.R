library(tidyverse)
library(ashr)


# Load the GWAS results ---------------------------------------------------
gwas <- list.files("data/gemma/output", "*\\.assoc\\.txt", full.names = TRUE) %>%
  map(read_tsv)

# Location-year-specific degrees of freedom
df <- list.files("data/gemma/output", "single[0-9]{1,2}\\.log\\.txt", full.names = TRUE) %>%
  map_int(function(f) {
    temp <- read_lines(f)
    temp[str_detect(temp, "analyzed individuals")] %>%
      str_split(., " ") %>%
      sapply(., function(x) x[7]) %>%
      as.integer()
  })
df <- df - 9 - 1

# Add location-year names to the data matrices
sites <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  pull(Site) %>%
  unique()

names(gwas) <- sites
names(df) <- sites


# Adaptive shrinkage ------------------------------------------------------
# Computation of posterior effect and standard error estimates
ash_res <- map2(gwas, df, function(g, d) {
  ash(g$beta, g$se, mixcompdist = "uniform", df = d, outputlevel = 2, 
      control = list(maxiter.sqp = 5000))
})

# Summarize significance by local false sign rate
sval_05 <- map_int(ash_res, function(a) as.integer(sum(get_svalue(a) <= 0.05)))
sval_01 <- map_int(ash_res, function(a) as.integer(sum(get_svalue(a) <= 0.01)))

tibble(Site = sites, Sval0.05 = sval_05, Sval0.01 = sval_01) %>%
  gather(Threshold, Sig, Sval0.05, Sval0.01) %>%
  mutate(Threshold = str_replace(Threshold, "Sval", "")) %>%
  ggplot(., aes(x = Site, y = Sig, group = Threshold)) + theme_classic() +
    geom_col(aes(fill = Threshold), colour = "black", position = position_dodge()) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(x = "", y = "# Significant SNPs")
ggsave("figures/single/sig_snps.pdf", width = 12, height = 6, units = "in", dpi = 300)

# Compute 95% credible intervals for SNPs at s-value <= 0.05
ash_ci <- map2(ash_res, sites, function(a, s) {
  temp <- ashci(a, lfsr_threshold = 0.05)
  rownames(temp) <- gwas[[s]]$rs
  temp[!is.na(temp[, 1]) & !is.na(temp[, 2]), ]
})

write_rds(list(ash = ash_res, ci = ash_ci, idx01 = sval_01, idx05 = sval_05), 
          "data/gemma/ash_results.rds")
