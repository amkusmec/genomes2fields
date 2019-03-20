library(tidyverse)
library(mgcv)


# Keep plot records only from sites with weather data ---------------------
yield <- read_rds("data/phenotype/yield_augment.rds")


# Identify sites for imputation -------------------------------------------
missing <- yield %>%
  group_by(Site) %>%
  summarise(Stand = sum(is.na(Stand))/n(), 
            StalkLodging = sum(is.na(StalkLodging))/n(), 
            RootLodging = sum(is.na(RootLodging))/n()) %>%
  ungroup()

# Identify sites without 100% coverage
stand_0 <- missing %>%
  filter(Stand == 0) %>%
  pull(Site)
stalk_0 <- missing %>%
  filter(StalkLodging == 0) %>%
  pull(Site)
root_0 <- missing %>%
  filter(RootLodging == 0) %>%
  pull(Site)

# Identify sites with at least 40% coverage
stand_40 <- missing %>%
  filter(Stand <= 0.6) %>%
  pull(Site)
stalk_40 <- missing %>%
  filter(StalkLodging <= 0.6) %>%
  pull(Site)
root_40 <- missing %>%
  filter(RootLodging <= 0.6) %>%
  pull(Site)


# Modified yield table for 100% coverage sites ----------------------------
yield %>%
  mutate(Stand = if_else(Site %in% stand_0, Stand, as.integer(NA)), 
         StalkLodging = if_else(Site %in% stalk_0, StalkLodging, as.integer(NA)), 
         RootLodging = if_else(Site %in% root_0, RootLodging, as.integer(NA))) %>%
  write_rds(., "data/phenotype/yield_agron0.rds")


# Imputation for 40% coverage sites ---------------------------------------
yield <- yield %>%
  mutate(Stand = if_else(Site %in% stand_40, Stand, as.integer(NA)), 
         StalkLodging = if_else(Site %in% stalk_40, StalkLodging, as.integer(NA)), 
         RootLodging = if_else(Site %in% root_40, RootLodging, as.integer(NA)))

# Remove sites without complete field layout information
no_field <- yield %>%
  group_by(Site) %>%
  summarise(Row = sum(is.na(Row))/n(), 
            Column = sum(is.na(Column))/n()) %>%
  ungroup() %>%
  filter(Row > 0, Column > 0) %>%
  pull(Site)

stand_40_n <- intersect(stand_40, setdiff(no_field, stand_0))
stalk_40_n <- intersect(stalk_40, setdiff(no_field, stalk_0))
root_40_n <- intersect(root_40, setdiff(no_field, root_0))

stand_40 <- setdiff(stand_40, union(stand_0, no_field))
stalk_40 <- setdiff(stalk_40, union(stalk_0, no_field))
root_40 <- setdiff(root_40, union(root_0, no_field))

yield <- yield %>%
  mutate(Stand = if_else(Site %in% stand_40_n, as.integer(NA), Stand), 
         StalkLodging = if_else(Site %in% stalk_40_n, as.integer(NA), StalkLodging), 
         RootLodging = if_else(Site %in% root_40_n, as.integer(NA), RootLodging))

# Impute stand count
for (i in stand_40) {
  cat(i, "\n")
  
  miss_idx <- which(is.na(yield$Stand) & yield$Site == i)
  if (length(miss_idx) == 0) next
  
  temp <- filter(yield, Site == i) %>%
    mutate(Replicate = as.factor(Replicate), 
           Block = as.factor(Block), 
           Rowf = as.factor(Row), 
           Colf = as.factor(Column))
  if (str_detect(i, "201[45]")) {
    temp_fit <- gam(Stand ~ Replicate + s(Replicate, Block, bs = "re") +
                      s(Rowf, bs = "re") + s(Colf, bs = "re") +
                      te(Row, Column, bs = "ps"), family = poisson(), 
                    data = filter(temp, !is.na(Stand)), method = "REML", 
                    optimizer = "outer", select = TRUE, 
                    drop.unused.levels = FALSE)
  } else if (str_detect(i, "201[67]")) {
    temp_fit <- gam(Stand ~ Replicate + s(Rowf, bs = "re") + 
                      s(Colf, bs = "re") + te(Row, Column, bs = "ps"), 
                    family = poisson(), data = filter(temp, !is.na(Stand)), 
                    method = "REML", optimizer = "outer", select = TRUE, 
                    drop.unused.levels = FALSE)
  }
  
  temp_pred <- predict.gam(temp_fit, filter(temp, is.na(Stand)), 
                           type = "response")
  yield$Stand[miss_idx] <- temp_pred
}

# Impute stalk lodging
for (i in stalk_40) {
  cat(i, "\n")
  
  miss_idx <- which(is.na(yield$StalkLodging) & yield$Site == i)
  if (length(miss_idx) == 0) next
  
  temp <- filter(yield, Site == i) %>%
    mutate(Replicate = as.factor(Replicate), 
           Block = as.factor(Block), 
           Rowf = as.factor(Row), 
           Colf = as.factor(Column))
  if (i == "GAH1_2015") {
    temp_fit <- gam(StalkLodging ~ s(Rowf, bs = "re") + s(Colf, bs = "re") +
                      te(Row, Column, bs = "ps"), family = poisson(), 
                    data = filter(temp, !is.na(StalkLodging)), method = "REML", 
                    optimizer = "outer", select = TRUE, drop.unused.levels = TRUE)
  } else if (str_detect(i, "201[45]")) {
    temp_fit <- gam(StalkLodging ~ Replicate + s(Replicate, Block, bs = "re") +
                      s(Rowf, bs = "re") + s(Colf, bs = "re") +
                      te(Row, Column, bs = "ps"), family = poisson(), 
                    data = filter(temp, !is.na(StalkLodging)), method = "REML", 
                    optimizer = "outer", select = TRUE, 
                    drop.unused.levels = FALSE)
  } else if (str_detect(i, "201[67]")) {
    temp_fit <- gam(StalkLodging ~ Replicate + s(Rowf, bs = "re") + 
                      s(Colf, bs = "re") + te(Row, Column, bs = "ps"), 
                    family = poisson(), data = filter(temp, !is.na(StalkLodging)), 
                    method = "REML", optimizer = "outer", select = TRUE, 
                    drop.unused.levels = FALSE)
  }
  
  temp_pred <- predict.gam(temp_fit, filter(temp, is.na(StalkLodging)), 
                           type = "response")
  yield$StalkLodging[miss_idx] <- temp_pred
}

# Impute root lodging
for (i in root_40) {
  cat(i, "\n")
  
  miss_idx <- which(is.na(yield$RootLodging) & yield$Site == i)
  if (length(miss_idx) == 0) next
  
  temp <- filter(yield, Site == i) %>%
    mutate(Replicate = as.factor(Replicate), 
           Block = as.factor(Block), 
           Rowf = as.factor(Row), 
           Colf = as.factor(Column))
  if (i == "GAH1_2015") {
    temp_fit <- gam(RootLodging ~ s(Rowf, bs = "re") + s(Colf, bs = "re") +
                      te(Row, Column, bs = "ps"), family = poisson(), 
                    data = filter(temp, !is.na(RootLodging)), method = "REML", 
                    optimizer = "outer", select = TRUE, drop.unused.levels = TRUE)
  } else if (str_detect(i, "201[45]")) {
    temp_fit <- gam(RootLodging ~ Replicate + s(Replicate, Block, bs = "re") +
                      s(Rowf, bs = "re") + s(Colf, bs = "re") +
                      te(Row, Column, bs = "ps"), family = poisson(), 
                    data = filter(temp, !is.na(RootLodging)), method = "REML", 
                    optimizer = "outer", select = TRUE, 
                    drop.unused.levels = FALSE)
  } else if (str_detect(i, "201[67]")) {
    temp_fit <- gam(RootLodging ~ Replicate + s(Rowf, bs = "re") + 
                      s(Colf, bs = "re") + te(Row, Column, bs = "ps"), 
                    family = poisson(), data = filter(temp, !is.na(RootLodging)), 
                    method = "REML", optimizer = "outer", select = TRUE, 
                    drop.unused.levels = FALSE)
  }
  
  temp_pred <- predict.gam(temp_fit, filter(temp, is.na(RootLodging)), 
                           type = "response")
  yield$RootLodging[miss_idx] <- temp_pred
}

write_rds(yield, "data/phenotype/yield_agron40.rds")
