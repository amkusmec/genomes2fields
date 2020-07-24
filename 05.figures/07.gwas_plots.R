library(tidyverse)
library(grid)
library(gridExtra)
source("src/manhattan_gemma.R")
source("src/qqplot_gemma.R")


gwas <- list.files("data/gemma/output", "*\\.assoc\\.txt", full.names = TRUE) %>%
  map(read_tsv)
phenotypes <- list.files("data/gemma/output", "*\\.assoc\\.txt") %>%
  str_remove(., "\\.assoc\\.txt") %>%
  str_remove(., "norm_")

m3 <- read_rds("data/gemma/norm_snp_mash_all.rds")

lay <- matrix(c(1, 1, 1, 2, 2), nrow = 1)

for (i in seq_along(phenotypes)) {
  pA <- gwas[[i]] %>%
    mutate(pval = m3$result$lfsr[, i]) %>%
    manhattan_gemma(., cutoff = -log10(0.1), ylab = expression(-log[10](s-value)))
  pB <- qqplot_gemma(gwas[[i]]$p_wald)
  gA <- grobTree(ggplotGrob(pA), 
                 textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                          hjust = "left", vjust = "top", 
                          gp = gpar(fontface = "bold", fontsize = 14)))
  gB <- grobTree(ggplotGrob(pB), 
                 textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                          hjust = "left", vjust = "top", 
                          gp = gpar(fontface = "bold", fontsize = 14)))
  gp <- arrangeGrob(gA, gB, layout_matrix = lay)
  ggsave(paste0("figures/select/gwas_plots/", phenotypes[i], ".png"), plot = gp, 
         width = 12, height = 4, units = "in", dpi = 300)
}
