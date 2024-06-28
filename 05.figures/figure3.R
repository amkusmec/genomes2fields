library(tidyverse)
library(grid)
library(gridExtra)

sel_vars <- read_rds("figures/ceris/criteria.rds")

sel_varsX <- sel_vars %>% 
  select(Method:Site, eBIC) %>% 
  group_by(Site) %>% 
  mutate(Best = if_else(eBIC == min(eBIC), TRUE, FALSE)) %>% 
  ungroup() %>% 
  mutate(Category = if_else(Site == "Full", "Full", 
                            if_else(Best & Site != "Full", "Site", "No")))
sel_varsX0 <- filter(sel_varsX, Site == "Full") %>% 
  mutate(X = 1:3 - 0.5, XEND = 1:3 + 0.5)
sel_varsX1 <- filter(sel_varsX, Site != "Full") %>% 
  mutate(Method2 = case_when(
    Method == "CERIS"   ~ 1, 
    Method == "climwin" ~ 2, 
    Method == "GA"      ~ 3
  ))

evar_lbls <- function(y) {
  sapply(y, function(x) {
    if (x == "NET") {
      expression(bold("NET (mm)"))
    } else if (x == "PPT") {
      expression(bold("PPT (mm)")) 
    } else if (x == "SR") {
      expression(bold("SR (W "*m^-2*")"))
    } else if (x == "TMAX") {
      expression(bold("TMAX ("*degree*"C)"))
    } else if (x == "TMIN") {
      expression(bold("TMIN ("*degree*"C)"))
    }
  })
}

pB <- ggplot() + theme_bw() + 
  geom_line(aes(x = Method2, y = eBIC, group = Site), sel_varsX1, 
            colour = "grey70", alpha = 0.5, 
            position = position_dodge(width = 0.25)) + 
  geom_segment(aes(x = X, xend = XEND, y = eBIC, yend = eBIC), sel_varsX0, 
               colour = "red", linewidth = 1) + 
  geom_point(aes(x = Method2, y = eBIC, group = Site, colour = Category), 
             sel_varsX1, position = position_dodge(width = 0.25), size = 1, 
             alpha = 0.5) + 
  scale_colour_manual(values = c("Site" = "green", "No" = "black"), 
                      labels = c("Site" = "Yes", "No" = "No")) + 
  scale_x_continuous(breaks = 1:3, labels = c("CERIS", "climwin", "GA")) + 
  labs(x = "", colour = "Best model?", tag = "(b)") + 
  theme(legend.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent"), 
        legend.position = "inside", 
        legend.position.inside = c(0.65, 0.8), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        legend.key.spacing.y = unit(2, "points"), 
        legend.key.height = unit(0.5, "lines"), 
        axis.title = element_text(size = 8), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1))


bars <- sel_vars %>% 
  select(Method:Variables) %>% 
  unnest(cols = Variables) %>% 
  separate("Variables", c("Variable", "Start", "End"), sep = "_") %>% 
  mutate(across(Start:End, ~ as.numeric(str_remove(.x, "X")))) %>% 
  mutate(Site = factor(Site, ordered = TRUE, 
                       levels = c(rev(sort(setdiff(unique(Site), "Full"))), "Full"))) %>% 
  group_by(Method, Site, Variable) %>% 
  mutate(Multi = as.character(1:n())) %>% 
  ungroup()

pA <- ggplot(bars) + theme_bw() + 
  geom_vline(xintercept = 1, linetype = 1, colour = "yellow2", linewidth = 1.2) + 
  geom_segment(aes(x = Start, xend = End, y = Site, yend = Site, group = Site), 
               linewidth = 0.25) + 
  geom_point(aes(x = Start, y = Site, colour = Multi, shape = Multi), size = 1) + 
  geom_point(aes(x = End, y = Site, colour = Multi, shape = Multi), size = 1) + 
  scale_colour_manual(values = c("1" = "black", "2" = "red", "3" = "blue", "4" = "green")) + 
  scale_shape_manual(values = c("1" = 19, "2" = 17, "3" = 15, "4" = 18)) + 
  scale_x_continuous(labels = scales::label_percent()) + 
  facet_grid(Method ~ Variable, labeller = labeller(Variable = as_labeller(evar_lbls, default = label_parsed))) + 
  guides(colour = "none", shape = "none") + 
  labs(x = "% CHU to anthesis", y = "Environment", tag = "(a)") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 5), 
        axis.title = element_text(size = 8), 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(face = "bold", size = 8), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(colour = "white"), 
        panel.grid.minor.x = element_line(colour = "white"), 
        panel.background = element_rect(fill = "grey95"))

blank_grob <- ggplot() + theme(panel.background = element_rect(fill = "white"))

gp <- arrangeGrob(ggplotGrob(pA), ggplotGrob(pB), ggplotGrob(blank_grob), 
                  layout_matrix = matrix(c(1, 1, 1, 1, 2, 
                                           1, 1, 1, 1, 3), 
                                         nrow = 2, byrow = TRUE))
ggsave("figures/final/Figure_3.pdf", gp, width = 180, height = 90, units = "mm", dpi = 600)
