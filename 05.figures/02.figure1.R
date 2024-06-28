library(tidyverse)
library(grid)
library(gridExtra)


# Keep only hybrids that will be used for GxE analysis --------------------
taxa <- rownames(read_rds("data/gbs/add_snps.rds")$GD)

blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- rep(names(blue), times = sapply(blue, function(x) nrow(x$blue)))
blue <- blue %>%
  map_df(function(x) x$blue) %>%
  mutate(Site = sites) %>%
  filter(PedigreeNew %in% taxa) %>%
  select(Site, PedigreeNew, BLUE) %>%
  filter(Site != "NEH3_2015")
boot <- read_rds("data/phenotype/vcomp_bootstrap.rds")
boot <- boot[names(boot) != "NEH3_2015"]


# Across-site correlations ------------------------------------------------
sites <- unique(blue$Site) %>% sort()
n_sites <- length(sites)
cor_mat <- matrix(0, nrow = n_sites, ncol = n_sites)
diag(cor_mat) <- NA
sig_mat <- cor_mat

for (i in 1:(n_sites - 1)) {
  for (j in (i + 1):n_sites) {
    # Above-diagonal = number of common hybrids
    temp <- blue %>%
      filter(Site == sites[i] | Site == sites[j]) %>%
      spread(Site, BLUE, fill = NA)
    temp <- temp[!is.na(temp[[2]]) & !is.na(temp[[3]]), ]
    cor_mat[i, j] <- nrow(temp)
    
    # Below-diagonal = Kendall's tau b
    cor_mat[j, i] <- cor(temp[[2]], temp[[3]], method = "kendall")
    
    # Significance test
    sig_mat[j, i] <- cor.test(temp[[2]], temp[[3]], method = "kendall")$p.value
  }
}

dimnames(cor_mat) <- dimnames(sig_mat) <- list(sites, sites)


# Plot the correlation matrix ---------------------------------------------
cor_mat_lower <- cor_mat
cor_mat_lower[upper.tri(cor_mat_lower, diag = FALSE)] <- NA
cor_mat_lower <- cor_mat_lower %>%
  as_tibble(rownames = "Site1") %>%
  gather(Site2, Tau, -Site1, na.rm = TRUE)

cor_mat_upper <- cor_mat
cor_mat_upper[lower.tri(cor_mat_upper, diag = FALSE)] <- NA
cor_mat_upper <- cor_mat_upper %>%
  as_tibble(rownames = "Site1") %>%
  gather(Site2, N_hyb, -Site1, na.rm = TRUE)

sig_mat[upper.tri(sig_mat, diag = FALSE)] <- NA
sig_mat <- as_tibble(sig_mat, rownames = "Site1") %>%
  gather(Site2, P, -Site1, na.rm = TRUE) %>%
  mutate(Q = p.adjust(P, method = "fdr"), 
         Sig = if_else(Q <= 0.05, "*", ""))

cor_mat_lower <- cor_mat_lower %>% 
  mutate(across(Site1:Site2, ~ str_replace(.x, "_20", "\\(") %>% paste0(")")))
cor_mat_upper <- cor_mat_upper %>% 
  mutate(across(Site1:Site2, ~ str_replace(.x, "_20", "\\(") %>% paste0(")")))
sig_mat <- sig_mat %>% 
  mutate(across(Site1:Site2, ~ str_replace(.x, "_20", "\\(") %>% paste0(")")))

pA <- ggplot() + theme_classic() + 
  geom_tile(aes(x = Site2, y = Site1, fill = Tau), cor_mat_lower) +
  # geom_text(aes(x = Site2, y = Site1, label = N_hyb), cor_mat_upper, size = 2) +
  geom_point(aes(x = Site2, y = Site1, colour = N_hyb), cor_mat_upper, size = 1.5, shape = 15) + 
  geom_text(aes(x = Site2, y = Site1, label = Sig), sig_mat, size = 2.5, nudge_x = -0.5) +
  scale_fill_distiller(type = "div", palette = "RdBu", direction = 1) +
  scale_colour_steps(low = "#f7fcf5", high = "#005a32", show.limits = TRUE, 
                     breaks = seq(0, 350, 50)) + 
  labs(x = "", y = "", fill = expression(tau[b]), 
       colour = "Hybrids", tag = "(b)") + coord_flip() + 
  scale_x_discrete(position = "top") + scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(hjust = -0.0625, angle = 90, size = 4.5), 
        axis.text.y = element_text(size = 5), 
        legend.position = "bottom", 
        legend.key.width = unit(0.05, "npc"), 
        legend.key.height = unit(0.1, "npc"), 
        legend.text = element_text(size = 6, angle = 45, vjust = 0.25), 
        legend.title = element_text(size = 8), 
        legend.title.position = "top", 
        plot.margin = unit(c(5.5, 2, 5.5, 2), "points"))


# Heterogeneity in scaled genetic variance
gvar <- boot %>%
  map_df(function(l) {
    totals <- rowSums(l$components)
    tibble(VarG = median(l$components[, "u:PedigreeNew"]/totals, na.rm = TRUE),
           Lower = quantile(l$components[, "u:PedigreeNew"]/totals, probs = 0.025, na.rm = TRUE),
           Upper = quantile(l$components[, "u:PedigreeNew"]/totals, probs = 0.975, na.rm = TRUE))
  }) %>%
  mutate(Site = names(boot)) %>%
  dplyr::select(Site, everything())
gvar <- boot %>%
  map_df(function(l) {
    totals <- rowSums(l$components)
    temp <- l$components[, "u:PedigreeNew"]/totals
    s <- sd(temp, na.rm = TRUE)/sqrt(sum(!is.na(temp)))
    tibble(VarG = mean(temp, na.rm = TRUE), 
           Lower = VarG - 3*s, 
           Upper = VarG + 3*s)
  }) %>%
  mutate(Site = names(boot), 
         Lower = if_else(Lower < 0, 0, Lower)) %>%
  dplyr::select(Site, everything())

# pB <- ggplot(gvar, aes(x = Site)) + theme_classic() +
#   geom_hline(aes(yintercept = Y), colour = "grey80", 
#              data = tibble(Y = c(0, 0.25, 0.5, 0.75))) + 
#   geom_point(aes(y = VarG)) +
#   geom_segment(aes(x = Site, xend = Site, y = Lower, yend = Upper)) + 
#   scale_y_continuous(labels = scales::percent, limits = c(0, 0.9)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         plot.margin = unit(c(7.5, 5.5, 5.5, 3.5), units = "points")) +
#   labs(x = "", y = expression(paste("% ", sigma[g]^2)), tag = "a") + 
#     theme(axis.title = element_text(size = 15), 
#           axis.text = element_text(size = 12))

gvar <- gvar %>% separate("Site", c("Experiment", "Year"), sep = "_", remove = FALSE)
pB <- ggplot(gvar, aes(x = Experiment)) + theme_classic() +
  geom_hline(aes(yintercept = Y), colour = "grey80", 
             data = tibble(Y = c(0, 0.25, 0.5, 0.75))) + 
  geom_point(aes(y = VarG, group = Year, colour = Year), position = position_dodge(0.85)) +
  geom_linerange(aes(ymin = Lower, ymax = Upper, 
                   group = Year, colour = Year), position = position_dodge(0.85)) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.9)) +
  scale_colour_brewer(type = "qual", palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.margin = unit(c(7.5, 5.5, 5.5, 3.5), units = "points")) +
  labs(x = "", y = expression(paste("% ", sigma[g]^2)), tag = "(a)", colour = "") + 
  theme(axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        legend.position = "inside", 
        legend.position.inside = c(0.45, 0.9), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8), 
        legend.background = element_rect(fill = "transparent"), 
        legend.direction = "horizontal")

lay <- matrix(c(1, 1, 2, 2, 2), nrow = 1)
# grob1 <- grobTree(ggplotGrob(pA), 
#                   textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
#                            hjust = "left", vjust = "top", 
#                            gp = gpar(fontface = "bold", fontsize = 14)))
# grob2 <- grobTree(ggplotGrob(pB), 
#                   textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
#                            hjust = "left", vjust = "top", 
#                            gp = gpar(fontface = "bold", fontsize = 14)))
# gp <- arrangeGrob(grob2, grob1, layout_matrix = lay)
gp <- arrangeGrob(ggplotGrob(pB), ggplotGrob(pA), layout_matrix = lay)
# ggsave("figures/single/gei_summary_miniscules.pdf", gp, width = 17, height = 8, units = "in", 
       # dpi = 300)
png("figures/final/figure01.png", width = 18, height = 9, 
    units = "in", res = 200)
plot(gp)
dev.off()

{
  gp <- arrangeGrob(ggplotGrob(pB), ggplotGrob(pA), layout_matrix = matrix(c(1, 1, 2, 2, 2), ncol = 1))
  ggsave("figures/final/Figure_1.pdf", gp, width = 80, height = 200, units = "mm", dpi = 600)
}
