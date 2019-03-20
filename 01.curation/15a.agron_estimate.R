library(tidyverse)
library(purrrlyr)
library(mgcv)


# Keep plot records only from sites with weather data ---------------------
yield <- read_rds("data/phenotype/yield_augment.rds")

cor(yield[, c("Yield", "Stand", "StalkLodging", "RootLodging")], use = "complete")


# Assess the pervasiveness of missing agronomic parameters ----------------
missing <- yield %>%
  group_by(Site) %>%
  summarise(Stand = sum(is.na(Stand))/n(), 
            StalkLodging = sum(is.na(StalkLodging))/n(), 
            RootLodging = sum(is.na(RootLodging))/n(), 
            n = n()) %>%
  ungroup()

missing %>%
  gather(Parameter, Value, Stand:RootLodging) %>%
  ggplot(., aes(x = Value, fill = Parameter)) + theme_bw() +
    geom_histogram(binwidth = 0.05, alpha = 0.8, colour = "black") +
    scale_fill_manual(values = c("Stand" = "red", "StalkLodging" = "green", 
                                 "RootLodging" = "blue")) +
    facet_wrap(~ Parameter) + guides(fill = FALSE) + 
    labs(x = "Proportion Missing Plots", y = "Count") +
    ggtitle("N = 64 location-years")
ggsave("figures/munge/missing_agronomy.pdf", width = 8, height = 5, units = "in", dpi = 300)


# An experiment to assess different amounts of missing data ---------------
temp <- yield %>%
  filter(Site == "NYH3_2015") %>%
  mutate(Rowf = factor(Row), 
         Colf = factor(Column), 
         Replicate = factor(Replicate), 
         Block = factor(Block),
         Pedigree = factor(Pedigree))

ggplot(temp, aes(x = Column, y = Row, fill = StalkLodging)) + theme_bw() +
  geom_tile() + labs(x = "Column", y = "Row", fill = "Stalk Lodging") +
  scale_fill_distiller(type = "seq", palette = "OrRd", direction = 1)
ggsave("figures/munge/stalk_NYH3_2015.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Missing at random -------------------------------------------------------
rand_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(nrow(temp), floor(r$Missing[1]*nrow(temp)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

rand_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 10, linetype = 2, colour = "red") +
    labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() + ggtitle("Missing at random")
ggsave("figures/munge/spatial_missing_random.pdf", width = 7, height = 5, units = "in", dpi = 300)


# Different quarters missing ----------------------------------------------
# NE corner
corner_idx <- with(temp, which(Row >= 9 &Column >= 15))
corner_ne_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(corner_idx, floor(r$Missing[1]*length(corner_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

corner_ne_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 10, linetype = 2, colour = "red") +
    labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() + ggtitle("Missing NE corner")
ggsave("figures/munge/spatial_missing_corner_ne.pdf", width = 7, height = 5, units = "in", dpi = 300)

# SE corner
corner_idx <- with(temp, which(Row >= 9 & Column <= 15))
corner_se_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(corner_idx, floor(r$Missing[1]*length(corner_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

corner_se_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 10, linetype = 2, colour = "red") +
    labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() + ggtitle("Missing SE corner")
ggsave("figures/munge/spatial_missing_corner_se.pdf", width = 7, height = 5, units = "in", dpi = 300)

# SW corner
corner_idx <- with(temp, which(Row <= 9 & Column <= 15))
corner_sw_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(corner_idx, floor(r$Missing[1]*length(corner_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

corner_sw_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 10, linetype = 2, colour = "red") +
    labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() + ggtitle("Missing SW corner")
ggsave("figures/munge/spatial_missing_corner_sw.pdf", width = 7, height = 5, units = "in", dpi = 300)

# NW corner
corner_idx <- with(temp, which(Row >= 9 & Column <= 15))
corner_nw_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(corner_idx, floor(r$Missing[1]*length(corner_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

corner_nw_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 10, linetype = 2, colour = "red") +
    labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() + ggtitle("Missing NW corner")
ggsave("figures/munge/spatial_missing_corner_nw.pdf", width = 7, height = 5, units = "in", dpi = 300)


# Different halves missing ------------------------------------------------
# Replicate 1 missing
rep_idx <- with(temp, which(Replicate == "1"))
rep1_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(rep_idx, floor(r$Missing[1]*length(rep_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

rep1_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_hline(yintercept = 10, linetype = 2, colour = "red") +
    labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() + ggtitle("Missing replicate 1")
ggsave("figures/munge/spatial_missing_rep1.pdf", width = 7, height = 5, units = "in", dpi = 300)

# Replicate 2 missing
rep_idx <- with(temp, which(Replicate == "2"))
rep2_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(rep_idx, floor(r$Missing[1]*length(rep_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

rep2_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Missing replicate 2")
ggsave("figures/munge/spatial_missing_rep2.pdf", width = 7, height = 5, units = "in", dpi = 300)

# E half missing
rep_idx <- with(temp, which(Column >= 16))
half_e_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(rep_idx, floor(r$Missing[1]*length(rep_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

half_e_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Missing E half")
ggsave("figures/munge/spatial_missing_half_e.pdf", width = 7, height = 5, units = "in", dpi = 300)

# W half missing
rep_idx <- with(temp, which(Column <= 15))
half_w_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(rep_idx, floor(r$Missing[1]*length(rep_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

half_w_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Missing W half")
ggsave("figures/munge/spatial_missing_half_w.pdf", width = 7, height = 5, units = "in", dpi = 300)


# Multiple corners missing ------------------------------------------------
# NW/SE corners missing
diag_idx <- with(temp, c(which(Row >= 9 & Column <= 15), 
                         which(Row <= 9 & Column >= 16)))
corner1_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(diag_idx, floor(r$Missing[1]*length(diag_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

corner1_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Missing NW/SE corners")
ggsave("figures/munge/spatial_missing_corner1.pdf", width = 7, height = 5, units = "in", dpi = 300)

# SW/NE corners missing
diag_idx <- with(temp, c(which(Row >= 9 & Column >= 16), 
                         which(Row <= 9 & Column <= 15)))
corner2_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(diag_idx, floor(r$Missing[1]*length(diag_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

corner2_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Missing NE/SW corners")
ggsave("figures/munge/spatial_missing_corner2.pdf", width = 7, height = 5, units = "in", dpi = 300)

# Only NE
diag_idx <- setdiff(1:nrow(temp), with(temp, which(Row >= 9 & Column >= 16)))
quarter1_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(diag_idx, floor(r$Missing[1]*length(diag_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

quarter1_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Only NE corner")
ggsave("figures/munge/spatial_missing_kept_ne.pdf", width = 7, height = 5, units = "in", dpi = 300)

# Only SE
diag_idx <- setdiff(1:nrow(temp), with(temp, which(Row <= 9 & Column >= 16)))
quarter2_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(diag_idx, floor(r$Missing[1]*length(diag_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

quarter2_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Only SE corner")
ggsave("figures/munge/spatial_missing_kept_se.pdf", width = 7, height = 5, units = "in", dpi = 300)

# Only SW
diag_idx <- setdiff(1:nrow(temp), with(temp, which(Row <= 9 & Column <= 15)))
quarter3_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(diag_idx, floor(r$Missing[1]*length(diag_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

quarter3_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Only SW corner")
ggsave("figures/munge/spatial_missing_kept_sw.pdf", width = 7, height = 5, units = "in", dpi = 300)

# Only NW
diag_idx <- setdiff(1:nrow(temp), with(temp, which(Row >= 9 & Column <= 15)))
quarter4_rpmse <- tibble(Missing = seq(0.05, 0.8, 0.05)) %>%
  by_row(function(r) {
    1:100 %>%
      map_dbl(function(x) {
        cat(r$Missing[1], x, "\n")
        idx <- sample(diag_idx, floor(r$Missing[1]*length(diag_idx)))
        temp$Mask <- temp$StalkLodging
        temp$Mask[idx] <- NA
        
        temp_fit <- gam(Mask ~ Replicate + s(Replicate, Block, bs = "re") +
                          s(Rowf, bs = "re") + s(Colf, bs = "re") + 
                          te(Row, Column, bs = "ps"), family = poisson(), 
                        data = filter(temp, !is.na(Mask)), method = "REML", 
                        optimizer = "outer", select = TRUE, 
                        drop.unused.levels = FALSE)
        temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Mask)), 
                                 type = "response")
        
        sqrt(mean((temp$StalkLodging[idx] - temp_pred)^2))
      })
  }, .to = "RPMSE")

quarter4_rpmse %>% unnest(RPMSE) %>%
  ggplot(., aes(x = Missing, y = RPMSE, group = Missing)) + theme_bw() +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = 2, colour = "red") +
  labs(x = "% Missing Plots", y = expression(log[10](RPMSE))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.8, 0.1)) +
  theme(panel.grid.minor = element_blank()) +
  annotation_logticks() + ggtitle("Only NW corner")
ggsave("figures/munge/spatial_missing_kept_nw.pdf", width = 7, height = 5, units = "in", dpi = 300)
