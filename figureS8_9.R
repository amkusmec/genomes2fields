library(tidyverse)
library(mashr)
library(grid)
library(gridExtra)


m2 <- read_rds("data/gemma/norm_snp_mash.rds")
phenotypes <-c("Intercept", "Whole season net evapotranspiration", "Mid-season solar radiation", 
               "Pre-anthesis max. temp.", "Post-anthesis max. temp.")
fills <- c("Intercept" = "grey60", 
           "Whole season net evapotranspiration" = "brown", 
           "Mid-season solar radiation" = "orange", 
           "Pre-anthesis max. temp." = "red", 
           "Post-anthesis max. temp." = "red")


# ED_1 plots --------------------------------------------------------------
# The first two panels recapitulate part of Figure 4
cmat1 <- cov2cor(m2$fitted_g$Ulist$ED_1)
rownames(cmat1) <- colnames(cmat1) <- phenotypes
cmat1[upper.tri(cmat1, diag = FALSE)] <- NA

pA1 <- as_tibble(cmat1, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = phenotypes, ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(phenotypes), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() + 
    geom_tile() + scale_x_discrete(position = "top") + 
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625)) + 
    labs(x = "", y = "", fill = "", subtitle = "relative frequency = 36.2%")

svd1 <- svd(cov2cor(m2$fitted_g$Ulist$ED_1))
pB1 <- tibble(Phenotype = phenotypes, 
              V = svd1$v[, 1]/svd1$v[, 1][which.max(abs(svd1$v[, 1]))]) %>%
  # mutate(Variable = str_remove(Phenotype, "_X[01]\\.[0-9]{2,3}_X[01]\\.[0-9]{1,3}")) %>%
  ggplot(., aes(x = Phenotype, y = V)) + theme_classic() + 
    geom_col(aes(fill = Phenotype), colour = "black") + 
    labs(x = "", y = "", fill = "", 
         subtitle = paste0("Eigenvector 1, (PVE = ", 
                           round(svd1$d[1]^2/sum(svd1$d^2), 2)*100, "%)")) + 
    scale_fill_manual(values = fills) + 
    scale_y_continuous(limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill = FALSE)

pC1 <- tibble(Phenotype = phenotypes, 
              V = svd1$v[, 2]/svd1$v[, 2][which.max(abs(svd1$v[, 2]))]) %>%
  # mutate(Variable = str_remove(Phenotype, "_X[01]\\.[0-9]{2,3}_X[01]\\.[0-9]{1,3}")) %>%
  ggplot(., aes(x = Phenotype, y = V)) + theme_classic() + 
    geom_col(aes(fill = Phenotype), colour = "black") + 
    labs(x = "", y = "", fill = "", 
         subtitle = paste0("Eigenvector 2, (PVE = ", 
                           round(svd1$d[2]^2/sum(svd1$d^2), 2)*100, "%)")) + 
    scale_fill_manual(values = fills) +
    scale_y_continuous(limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill = FALSE)

gA1 <- grobTree(ggplotGrob(pA1), 
                textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                         hjust = "left", vjust = "top", 
                         gp = gpar(fontface = "bold", fontsize = 14)))
gB1 <- grobTree(ggplotGrob(pB1), 
                textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                         hjust = "left", vjust = "top", 
                         gp = gpar(fontface = "bold", fontsize = 14)))
gC1 <- grobTree(ggplotGrob(pC1), 
                textGrob("C", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                         hjust = "left", vjust = "top", 
                         gp = gpar(fontface = "bold", fontsize = 14)))

lay1 <- matrix(c(1, 1, 2, 3), nrow = 1)
gp1 <- arrangeGrob(gA1, gB1, gC1, layout_matrix = lay1)
ggsave("figures/select/mash_ED1_all.pdf", plot = gp1, width = 12, height = 4, 
       units = "in", dpi = 300)


# ED_2 plots --------------------------------------------------------------
cmat2 <- cov2cor(m2$fitted_g$Ulist$ED_2)
rownames(cmat2) <- colnames(cmat2) <- phenotypes
cmat2[upper.tri(cmat2, diag = FALSE)] <- NA

pA2 <- as_tibble(cmat2, rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R)) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = phenotypes, ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(phenotypes), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R)) + theme_classic() + 
    geom_tile() + scale_x_discrete(position = "top") + 
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625)) + 
    labs(x = "", y = "", fill = "", subtitle = "relative frequency = 11.7%")

svd2 <- svd(cov2cor(m2$fitted_g$Ulist$ED_2))
pB2 <- tibble(Phenotype = phenotypes, 
              V = svd2$v[, 1]/svd2$v[, 1][which.max(abs(svd2$v[, 1]))]) %>%
  # mutate(Variable = str_remove(Phenotype, "_X[01]\\.[0-9]{2,3}_X[01]\\.[0-9]{1,3}")) %>%
  ggplot(., aes(x = Phenotype, y = V)) + theme_classic() + 
    geom_col(aes(fill = Phenotype), colour = "black") + 
    labs(x = "", y = "", fill = "", 
         subtitle = paste0("Eigenvector 1, (PVE = ", 
                           round(svd2$d[1]^2/sum(svd2$d^2), 2)*100, "%)")) + 
    scale_fill_manual(values = fills) + 
    scale_y_continuous(limits = c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill = FALSE)

gA2 <- grobTree(ggplotGrob(pA2), 
                textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                         hjust = "left", vjust = "top", 
                         gp = gpar(fontface = "bold", fontsize = 14)))
gB2 <- grobTree(ggplotGrob(pB2), 
                textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                         hjust = "left", vjust = "top", 
                         gp = gpar(fontface = "bold", fontsize = 14)))

lay2 <- matrix(c(1, 1, 2, 2), nrow = 1)
gp2 <- arrangeGrob(gA2, gB2, layout_matrix = lay2)
ggsave("figures/select/mash_ED2_all.pdf", plot = gp2, width = 10, height = 4, 
       units = "in", dpi = 300)

