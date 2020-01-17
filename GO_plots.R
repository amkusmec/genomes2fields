library(tidyverse)
library(GOSemSim)


cl <- read_csv("data/gemma/closest_vs_ld_perm.csv") %>%
  filter(perm <= 0.05) %>%
  mutate(Background = "LD")
cg <- read_csv("data/gemma/closest_vs_genome_perm.csv") %>%
  filter(perm <= 0.05) %>%
  mutate(Background = "Genome")

zmGO <- godata("org.Zmays.eg.db", keytype = "GID", ont = "BP", computeIC = TRUE)
ic_sub <- zmGO@IC[names(zmGO@IC) %in% unique(c(cg$GO, cl$GO))] %>%
  enframe() %>%
  dplyr::rename(GO = name, IC = value)

common_terms <- intersect(cl$GO, cg$GO)
shared <- bind_rows(cl, cg) %>%
  filter(GO %in% common_terms) %>%
  inner_join(., ic_sub, by = "GO") %>%
  arrange(desc(IC)) %>%
  mutate(Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE))

pA <- shared %>%
  dplyr::slice(1:50) %>%
  mutate(perm = -log10(perm)) %>%
  ggplot(.) + theme_bw() + 
    geom_point(aes(y = Name, x = perm, size = IC, fill = Background), 
               shape = 21, colour = "black", stroke = 0.5) + 
    scale_fill_manual(values = c("Genome" = "blue", "LD" = "red")) +
    labs(y = "", x = expression(paste(-log[10], "(p-value)")), 
         fill = "Background", size = "Information\nContent") + 
    guides(size = guide_legend(override.aes = list(fill = "black")), 
           fill = guide_legend(override.aes = list(size = 3)))

pS1 <- shared %>%
  dplyr::slice(201:250) %>%
  mutate(perm = -log10(perm)) %>%
  ggplot(.) + theme_bw() + 
    geom_point(aes(y = Name, x = perm, size = IC, fill = Background), 
               shape = 21, colour = "black", stroke = 0.5) + 
    scale_fill_manual(values = c("Genome" = "blue", "LD" = "red")) +
    labs(y = "", x = expression(paste(-log[10], "(p-value)")), 
         fill = "Background", size = "Information\nContent") + 
    guides(size = guide_legend(override.aes = list(fill = "black")), 
           fill = guide_legend(override.aes = list(size = 3)))

cl_unique <- cl %>%
  filter(!(GO %in% common_terms)) %>%
  inner_join(., ic_sub, by = "GO") %>%
  arrange(desc(IC), perm) %>%
  mutate(Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE))
cg_unique <- cg %>%
  filter(!(GO %in% common_terms)) %>%
  filter(!(GO %in% common_terms)) %>%
  inner_join(., ic_sub, by = "GO") %>%
  arrange(desc(IC), perm) %>%
  mutate(Name = factor(Name, levels = rev(unique(Name)), ordered = TRUE))

pS2a <- cl_unique %>%
  mutate(perm = -log10(perm)) %>%
  ggplot(.) + theme_bw() + 
    geom_point(aes(y = Name, x = perm, size = IC), shape = 21, colour = "black", 
               fill = "red", stroke = 0.5) + 
    labs(y = "", x = expression(paste(-log[10], "(p-value)")), 
         size = "Information\nContent")

pS2b <- cg_unique %>%
  mutate(perm = -log10(perm)) %>%
  ggplot(.) + theme_bw() + 
    geom_point(aes(y = Name, x = perm, size = IC), shape = 21, colour = "black", 
               fill = "blue", stroke = 0.5) + 
    labs(y = "", x = expression(paste(-log[10], "(p-value)")), 
         size = "Information\nContent")

coexpress <- read_rds("data/gemma/coexpress_enrich.rds")
mrna_go <- purrr::reduce(coexpress$Unique, union)
common <- union(intersect(mrna_go, cl$Name), intersect(mrna_go, cg$Name))
sapply(coexpress$Unique, function(l) sum(common %in% l))

library(VennDiagram)
term_venn <- draw.pairwise.venn(area1 = nrow(cl), area2 = nrow(cg), 
                                cross.area = length(common_terms), euler.d = FALSE, 
                                category = c("LD", "Genome"), scaled = FALSE, 
                                fill = c("red", "blue"), alpha = c(0.3, 0.3), 
                                cex = rep(2, 3), fontfamily = rep("sanserif", 3), 
                                cat.pos = rep(0, 2), cat.cex = rep(2, 2), 
                                cat.fontfamily = rep("sanserif", 2), 
                                cat.dist = rep(0.05, 2))


library(grid)
library(gridExtra)
gA <- grobTree(term_venn, 
               textGrob("A", x = unit(0.03, "npc"), y = unit(0.95, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gB <- grobTree(ggplotGrob(pA), 
               textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gp <- arrangeGrob(gA, gB, layout_matrix = matrix(c(1, 2, 2, NA, 2, 2), nrow = 2, byrow = TRUE))
ggsave("figures/select/go_main_figure.png", plot = gp, width = 12, height = 6, 
       units = "in", dpi = 300)

ggsave("figures/select/go_low_ic.pdf", plot = pS1, width = 8, height = 6, 
       units = "in", dpi = 300)

gS2a <- grobTree(ggplotGrob(pS2a), 
                 textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                          hjust = "left", vjust = "top", 
                          gp = gpar(fontface = "bold", fontsize = 14)))
gS2b <- grobTree(ggplotGrob(pS2b), 
                 textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                          hjust = "left", vjust = "top", 
                          gp = gpar(fontface = "bold", fontsize = 14)))
gp <- arrangeGrob(gS2a, gS2b, layout_matrix = matrix(c(1, 2), nrow = 1))
ggsave("figures/select/go_unique.pdf", plot = gp, width = 14, height = 6, 
       units = "in", dpi = 300)
