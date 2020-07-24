library(tidyverse)
library(grid)
library(gridExtra)


pc <- read_rds("data/gbs/pca.rds")

pA <- enframe(pc$sdev) %>%
  mutate(P = value^2/sum(value^2)) %>%
  slice(1:20) %>%
  ggplot(., aes(x = name, y = P)) + theme_bw() + 
    geom_line(linetype = 2) + geom_point() + 
    geom_hline(yintercept = 0.02, colour = "red", linetype = 2) + 
    labs(x = "PC", y = "% Variance") + 
    scale_y_continuous(labels = scales::percent)

scores <- as_tibble(pc$x, rownames = "PedigreeNew")

pB <- ggplot(scores, aes(x = PC1, y = PC2)) + theme_bw() + 
  geom_point(alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = 2, colour = "red") + 
  geom_vline(xintercept = 0, linetype = 2, colour = "red") + 
  labs(x = "PC1 (16.7%)", y = "PC2 (5.8%)")

pC <- ggplot(scores, aes(x = PC2, y = PC4)) + theme_bw() + 
  geom_point(alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = 2, colour = "red") + 
  geom_vline(xintercept = 0, linetype = 2, colour = "red") + 
  labs(x = "PC2 (5.8%)", y = "PC4 (4.3%)")

pD <- ggplot(scores, aes(x = PC3, y = PC4)) + theme_bw() + 
  geom_point(alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = 2, colour = "red") + 
  geom_vline(xintercept = 0, linetype = 2, colour = "red") + 
  labs(x = "PC3 (5.1%)", y = "PC4 (4.3%)")


gA <- grobTree(ggplotGrob(pA), 
               textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gB <- grobTree(ggplotGrob(pB), 
               textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gC <- grobTree(ggplotGrob(pC), 
               textGrob("C", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gD <- grobTree(ggplotGrob(pD), 
               textGrob("D", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))

gp <- arrangeGrob(gA, gB, gC, gD, layout_matrix = matrix(1:4, nrow = 2, byrow = TRUE))
ggsave("figures/munge/pca_plots.pdf", plot = gp, width = 10, height = 8, 
       units = "in", dpi = 300)
