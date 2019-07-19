library(tidyverse)


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
    tibble(Site = sites[i], Samples = txt[1], PVE = txt[3], SE = txt[4])
  }) %>%
  mutate_at(c("Samples"), as.integer) %>%
  mutate_at(c("PVE", "SE"), as.numeric)


# Plot sample size per site -----------------------------------------------
info %>% 
  select(Site, Samples) %>%
  rename(`Sample Size` = Samples) %>%
  ggplot(., aes(x = Site, y = `Sample Size`)) + theme_bw() +
    geom_point() + labs(x = "", y = "Sample Size") +
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/single/gwas_plots/common_samples.pdf", width = 8, height = 5, 
       units = "in", dpi = 300)


# Plot PVE for the relationship matrices ----------------------------------
info %>%
  mutate(Upper = PVE + SE,
         Upper = if_else(Upper > 1, 1, Upper), 
         Lower = PVE - SE, 
         Lower = if_else(Lower < 0, 0, Lower), 
         Shape = if_else(near(0, PVE, tol = 1e-4), "2", "1")) %>%
  ggplot(., aes(x = Site, y = PVE)) + theme_bw() +
    geom_pointrange(aes(ymin = Lower, ymax = Upper, shape = Shape), size = 1) +
    labs(x = "", y = "Percent Variance Explained") + guides(shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1))
ggsave("figures/single/gwas_plots/common_pve.pdf", width = 10, height = 6, units = "in", dpi = 300)


# Relationships between sample size and PVE -------------------------------
ggplot(info, aes(x = Samples, y = PVE)) + theme_classic() + geom_point()
ggplot(info, aes(x = Samples, y = SE)) + theme_classic() + geom_point()
ggplot(info, aes(x = PVE, y = SE)) + theme_classic() + geom_point()
