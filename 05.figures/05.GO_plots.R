library(tidyverse)
library(GOSemSim)
library(grid)
library(gridExtra)
library(VennDiagram)

anno_sc <- c("#66c2a5", "#fc8d62", "#8da0cb")


# Load the enriched GO terms
cl <- read_csv("data/gemma/closest_vs_ld_perm.csv") %>%
  filter(perm <= 0.05) %>%
  mutate(Background = "CvL")
cg <- read_csv("data/gemma/closest_vs_genome_perm.csv") %>%
  filter(perm <= 0.05) %>%
  mutate(Background = "CvG")
lg <- read_csv("data/gemma/ld_vs_genome_perm.csv") %>% 
  filter(perm <= 0.05) %>% 
  mutate(Background = "LvG")

# Calculate information content scores
zmGO <- godata("org.Zmays.eg.db", keytype = "GID", ont = "BP", computeIC = TRUE)
ic_sub <- zmGO@IC[names(zmGO@IC) %in% unique(c(cg$GO, cl$GO, lg$GO))] %>%
  enframe() %>%
  dplyr::rename(GO = name, IC = value)


### Various partitions
common3 <- intersect(cl$GO, intersect(cg$GO, lg$GO))
shared3 <- bind_rows(cl, cg, lg) %>%
  filter(GO %in% common3) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC)) %>%
  mutate(#Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE), 
         Combo = paste(unique(Background), collapse = "\n"))

common2_cl_cg <- setdiff(intersect(cl$GO, cg$GO), lg$GO)
shared2_cl_cg <- bind_rows(cl, cg) %>%
  filter(GO %in% common2_cl_cg) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC)) %>%
  mutate(#Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE), 
         Combo = paste(unique(Background), collapse = "\n"))

common2_cl_lg <- setdiff(intersect(cl$GO, lg$GO), cg$GO)
shared2_cl_lg <- bind_rows(cl, lg) %>%
  filter(GO %in% common2_cl_lg) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC)) %>%
  mutate(#Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE), 
         Combo = paste(unique(Background), collapse = "\n"))

common2_cg_lg <- setdiff(intersect(cg$GO, lg$GO), cl$GO)
shared2_cg_lg <- bind_rows(cg, lg) %>%
  filter(GO %in% common2_cg_lg) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC)) %>%
  mutate(#Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE), 
         Combo = paste(unique(Background), collapse = "\n"))


### Three-way venn diagram
term_venn <- draw.triple.venn(area1 = nrow(cl), area2 = nrow(cg), area3 = nrow(lg), 
                              n12 = length(intersect(cl$GO, cg$GO)), 
                              n23 = length(intersect(cg$GO, lg$GO)), 
                              n13 = length(intersect(cl$GO, lg$GO)), 
                              n123 = length(common3), 
                              euler.d = FALSE, 
                              category = c("Closest vs. LD", "Closest vs. genome", 
                                           "LD vs. genome"), scaled = FALSE, 
                              fill = anno_sc, alpha = rep(0.3, 3), 
                              cex = rep(1, 7), fontfamily = rep("sans", 7), 
                              cat.pos = c(345, 15, 180), cat.cex = rep(1, 3), 
                              cat.fontfamily = rep("sans", 3),
                              cat.fontface = rep("bold", 3), 
                              cat.dist = rep(0.05, 3), ind = FALSE)

# Common terms
df_a <- bind_rows(shared3, shared2_cg_lg, shared2_cl_lg, 
                  slice(shared2_cl_cg, 1:14)) %>%
  mutate(perm = -log10(perm), 
         IC = if_else(is.infinite(IC), max(IC[!is.infinite(IC)]), IC))
pA <- ggplot(df_a) + theme_bw() + 
  geom_point(aes(y = Name, x = IC, colour = perm, shape = Background), 
             position = position_jitter(height = 0.65, width = 0.35)) + 
  scale_colour_viridis_c() + 
    facet_wrap(~ Combo, ncol = 1, strip.position = "right", scales = "free_y") + 
    labs(y = "", x = "Information content", 
         fill = "Background", size = expression(paste(-log[10], "(p-value)")), tag = "(b)") + 
    # guides(size = guide_legend(override.aes = list(fill = "black")), 
    #        fill = guide_legend(override.aes = list(size = 3))) + 
  guides(shape = "none", colour = "none") + 
  theme(strip.text.y.right = element_text(size = 8, face = "bold", angle = 0), 
        strip.background = element_rect(fill = "white"), 
        axis.title = element_text(size = 8), 
        axis.text = element_text(size = 6))

gA <- grobTree(rectGrob(gp = gpar(fill = "white", col = "transparent")), 
               term_venn, 
               textGrob("(a)", x = unit(0.05, "npc"), y = unit(0.9, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontsize = 14)))


fake_shape <- tibble(Background = c("CvG", "CvL", "LvG"), 
                     X = 1, Y = 3:1)
fake_fill <- tibble(perm = seq(min(df_a$perm), max(df_a$perm), length.out = 20), 
                    X = 1:20, Y = 1:20)
blank_grob <- ggplot() + theme(panel.background = element_rect(fill = "white"))
leg_grob <- blank_grob + 
  geom_point(aes(x = X, y = Y, shape = Background), fake_shape, alpha = 0, size = 3) + 
  geom_point(aes(x = X, y = Y, colour = perm), fake_fill, alpha = 0, size = 3) + 
  scale_shape_manual(values = c("CvG" = 16, "CvL" = 17, "LvG" = 15)) + 
  scale_colour_viridis_c() +
  guides(shape = guide_legend(override.aes = list(alpha = 1)), 
         colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(colour = expression(paste(-log[10], "(p-value)"))) + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.box = "horizontal", 
        legend.background = element_rect(colour = "black"), 
        legend.position = "inside", 
        legend.position.inside = c(0.5, 0.5), 
        legend.text = element_text(size = 12), 
        legend.byrow = FALSE)

gp <- arrangeGrob(gA, ggplotGrob(pA), ggplotGrob(leg_grob), 
                  layout_matrix = matrix(c(1, 1, 3, 3, 2, 2, 2, 2), nrow = 2, byrow = TRUE))
ggsave("figures/final/Figure_6.pdf", gp, width = 180, height = 150, units = "mm", dpi = 600)



# pS1 <- shared %>%
#   dplyr::slice(201:250) %>%
#   mutate(perm = -log10(perm)) %>%
#   ggplot(.) + theme_bw() + 
#     geom_point(aes(y = Name, x = perm, size = IC, fill = Background), 
#                shape = 21, colour = "black", stroke = 0.5) + 
#     scale_fill_manual(values = c("Genome" = "blue", "LD" = "red")) +
#     labs(y = "", x = expression(paste(-log[10], "(p-value)")), 
#          fill = "Background", size = "Information\nContent") + 
#     guides(size = guide_legend(override.aes = list(fill = "black")), 
#            fill = guide_legend(override.aes = list(size = 3)))
# ggsave("figures/select/go_low_ic.pdf", plot = pS1, width = 8, height = 6, 
#        units = "in", dpi = 300)

# Uniquely enriched terms
cl_unique <- cl %>%
  filter(!(GO %in% c(common3, common2_cl_cg, common2_cl_lg))) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC), perm)
cg_unique <- cg %>%
  filter(!(GO %in% c(common3, common2_cl_cg, common2_cg_lg))) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC), perm)
lg_unique <- lg %>%
  filter(!(GO %in% c(common3, common2_cl_lg, common2_cg_lg))) %>%
  inner_join(ic_sub, by = "GO") %>%
  arrange(desc(IC), perm) %>% 
  slice(1:25)

pS2a <- cl_unique %>%
  mutate(perm = -log10(perm)) %>%
  ggplot() + theme_bw() + 
    geom_point(aes(y = Name, x = IC, size = perm), shape = 21, colour = "black", 
               fill = anno_sc[1], stroke = 0.5) + 
    labs(y = "", x = "Information content", tag = "A", 
         size = expression(paste(-log[10], "(p-value)")), 
         subtitle = "Closest vs. LD")

pS2b <- cg_unique %>%
  mutate(perm = -log10(perm)) %>%
  ggplot() + theme_bw() + 
  geom_point(aes(y = Name, x = IC, size = perm), shape = 21, colour = "black", 
             fill = anno_sc[2], stroke = 0.5) + 
  labs(y = "", x = "Information content", tag = "B", 
       size = expression(paste(-log[10], "(p-value)")), 
       subtitle = "Closest vs. genome")

pS2c <- lg_unique %>%
  mutate(perm = -log10(perm)) %>%
  ggplot() + theme_bw() + 
  geom_point(aes(y = Name, x = IC, size = perm), shape = 21, colour = "black", 
             fill = anno_sc[3], stroke = 0.5) + 
  labs(y = "", x = "Information content", tag = "C", 
       size = expression(paste(-log[10], "(p-value)")), 
       subtitle = "LD vs. genome")

gp <- arrangeGrob(ggplotGrob(pS2a), ggplotGrob(pS2b), ggplotGrob(pS2c), 
                  rectGrob(gp = gpar(fill = "white", col = "transparent")), 
                  layout_matrix = matrix(1:4, nrow = 2, byrow = TRUE))
png("figures/select/go_unique.png", width = 15, height = 9, units = "in", res = 300)
plot(gp)
dev.off()
