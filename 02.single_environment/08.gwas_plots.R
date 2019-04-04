library(tidyverse)
library(grid)
library(gridExtra)

source("src/manhattan_gemma.R")
source("src/qqplot_gemma.R")


# Manhattan and QQ plots --------------------------------------------------
blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- names(blue)
sites <- sites[!str_detect(sites, "2017$")]
rm(blue); gc()

for (i in seq_along(sites)) {
  assoc <- read_tsv(paste0("data/gemma/output/single", i, ".assoc.txt")) %>%
    rename(pval = p_lrt)
  
  g1 <- manhattan_gemma(assoc, cutoff = -log10(0.05/nrow(assoc)), title = sites[i])
  g2 <- qqplot_gemma(assoc$pval)
  
  g <- arrangeGrob(ggplotGrob(g1), ggplotGrob(g2), 
                   layout_matrix = matrix(c(1, 1, 2), nrow = 1))
  ggsave(paste0("figures/single/gwas_plots/single_", sites[i], ".pdf"), 
         g, width = 10, height = 5, units = "in", dpi = 300)
}
rm(assoc, g, g1, g2, i); gc()


# Scrape information from the log files -----------------------------------
# Collecting information on sample size, number of markers, percent variance
# explained by G, and standard error of PVE
info <- seq_along(sites) %>%
  map_df(function(i) {
    txt <- read_lines(paste0("data/gemma/output/single", i, ".log.txt"))
    txt <- txt[str_detect(txt, "analyzed") | str_detect(txt, "pve")]
    txt <- str_replace(txt, "## [:print:]* = ", "")
    tibble(Site = sites[i], Samples = txt[1], SNPs = txt[2], 
           PVE = txt[3], SE = txt[4])
  }) %>%
  mutate_at(c("Samples", "SNPs"), as.integer) %>%
  mutate_at(c("PVE", "SE"), as.numeric)


# Plot sample sizes and marker numbers ------------------------------------
info %>%
  select(Site:SNPs) %>%
  rename(`Sample Size` = Samples, `# SNPs` = SNPs) %>%
  gather(Variable, Value, -Site) %>%
  ggplot(., aes(x = Site, y = Value)) + theme_bw() +
    geom_point() + labs(x = "", y = "") +
    facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/single/gwas_plots/sample_snps.pdf", width = 8, height = 6, 
       units = "in", dpi = 300)


# Plot variance explained by the relationship matrix ----------------------
info %>%
  mutate(SE = if_else(near(0, PVE, tol = 1e-4), 0, SE), 
         Shape = if_else(near(0, PVE, tol = 1e-4), "2", "1")) %>%
  ggplot(., aes(x = Site, y = PVE)) + theme_bw() +
    geom_pointrange(aes(ymin = PVE - SE, ymax = PVE + SE, shape = Shape), size = 1) +
    labs(x = "", y = "Percent Variance Explained") + guides(shape = "none") +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    scale_y_continuous(labels = scales::percent, limits = c(-0.05, 1))
ggsave("figures/single/gwas_plots/pve.pdf", width = 10, height = 6, units = "in", dpi = 300)
