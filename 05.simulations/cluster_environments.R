library(tidyverse)
library(lubridate)
library(purrrlyr)

weather <- read_rds("data/weather/env_variables.rds") %>%
  mutate(Site = paste(Environment, Year, sep = "_"), 
         Date = paste(Year, Month, Day, sep = "-") %>% ymd())

yield <- read_rds("data/phenotype/yield_agron0.rds") %>%
  filter(Year < 2017) %>%
  select(Site, PedigreeNew, CHUA, Planted, Harvested) %>%
  group_by(Site) %>%
  summarise(CHUA = median(CHUA), 
            Planted = min(Planted), 
            Harvested = max(Harvested)) %>%
  ungroup()

thresholds <- seq(0.025, 1.5, 0.025)
dates <- yield %>%
  by_row(function(r) {
    temp <- weather %>%
      filter(Site == r$Site[1], Date >= r$Planted[1], Date <= r$Harvested)
    chu <- cumsum(temp$CHU)/r$CHUA
    indices <- map_int(thresholds, function(x) {
      as.integer(min(which(chu - x >= 0)))
    })
    ymd(temp$Date[indices])
  }, .collate = "cols", .labels = FALSE) %>%
  mutate_all(funs(as_date))
names(dates) <- make.names(thresholds)

dates <- dates %>%
  mutate(X0 = as.character(yield$Planted)) %>%
  select(X0, everything())

date_idx <- matrix(0L, nrow = nrow(dates), ncol = ncol(dates))
for (i in 1:ncol(date_idx)) {
  for (j in 1:nrow(date_idx)) {
    date_idx[j, i] <- which(weather$Site == yield$Site[j] &
                              weather$Date == dates[[i]][j])
  }
}
colnames(date_idx) <- names(dates)
rownames(date_idx) <- yield$Site

variables <- tibble(Var = c("NET", "NET", "TMAX", "TMIN", "TMIN"), 
                    Start = c("X0.025", "X0.65", "X0.775", "X0.05", "X0.15"), 
                    End = c("X0.45", "X1.275", "X0.875", "X1.5", "X0.65"))
yield <- yield %>%
  by_row(function(r) {
    c(sum(weather$NET[date_idx[r$Site, "X0.025"]:date_idx[r$Site, "X0.45"]]), 
      sum(weather$NET[date_idx[r$Site, "X0.65"]:date_idx[r$Site, "X1.275"]]), 
      mean(weather$TMAX[date_idx[r$Site, "X0.775"]:date_idx[r$Site, "X0.875"]]), 
      mean(weather$TMIN[date_idx[r$Site, "X0.05"]:date_idx[r$Site, "X1.5"]]), 
      mean(weather$TMIN[date_idx[r$Site, "X0.15"]:date_idx[r$Site, "X0.65"]]))
  }, .collate = "cols")
names(yield)[5:9] <- c("NET_X0.025_X0.45", "NET_X0.65_X1.275", "TMAX_X0.775_X0.875", 
                       "TMIN_X0.05_X1.5", "TMIN_X0.15_X0.65")

vars <- as.matrix(yield[, 5:9])
rownames(vars) <- yield$Site
vars <- scale(vars, center = TRUE, scale = TRUE)
k <- lapply(2:15, function(kk) kmeans(vars, centers = kk, 
                                      iter.max = 50, nstart = 50))


meta <- read_csv("data/metadata/metadata_clean.csv") %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% yield$Site) %>%
  arrange(Site)
meta$Cluster <- k[[3]]$cluster

states <- map_data("state")
ggplot(data = states) +
  geom_polygon(aes(x = long, y = lat, group = group), colour = "black", fill = "grey") +
  geom_point(aes(x = Longitude, y = Latitude, colour = factor(Year)), meta) +
  facet_wrap(~ Cluster) + coord_fixed(1.3) + labs(x = "", y = "", colour = "") +
  theme_bw()

as_tibble(vars, rownames = "Site") %>%
  mutate(Cluster = k[[3]]$cluster) %>%
  gather(Variable, Value, -Site, -Cluster) %>%
  ggplot(., aes(x = Variable, y = Site, fill = Value)) + theme_bw() +
    geom_tile() + scale_fill_distiller(type = "div", palette = "RdBu") +
    facet_wrap(~ Cluster, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

as_tibble(vars, rownames = "Site") %>%
  mutate(Cluster = k[[3]]$cluster) %>%
  split(., .$Cluster) %>%
  map_chr(function(df) {
    v <- as.matrix(df[, c(-1, -7)])
    d <- apply(v, 1, function(x) sqrt(sum((x - k[[3]]$centers[df$Cluster[1], ])^2)))
    df$Site[which.min(d)]
  })
