### Identify days that require imputation of minimum/maximum temperature and/or 
### precipitation. Daily values are retrieved from the PRISM project and assigned 
### to the nearest 4 km x 4 km grid based on great circle distance.

library(tidyverse)
library(purrrlyr)
library(lubridate)
library(prism)
library(parallel)


# Keep sites for which we have complete information -----------------------
weather <- read_rds("data/weather/weather_munged.rds") %>%
  select(-DEW, -SR)
yield <- read_rds("data/phenotype/yield_munged.rds")
meta <- read_csv("data/metadata/metadata_clean.csv")

weather_sites <- with(weather, paste(Environment, Year, sep = "_")) %>% unique()
yield_sites <- with(yield, paste(Environment, Year, sep = "_")) %>% unique()
meta_sites <- with(meta, paste(Environment, Year, sep = "_")) %>% unique()
sites_ym <- intersect(yield_sites, meta_sites)
sites <- intersect(weather_sites, sites_ym)

weather <- weather %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites)
yield <- yield %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites_ym)
meta <- meta %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites_ym)

# Some additional filtering based on the availability of planting/harvest dates
ph_sites <- yield %>%
  group_by(Year, Environment) %>%
  summarise(Planted = min(Planted, na.rm = TRUE), 
            Harvested = max(Harvested, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(as.character(Planted)), !is.na(as.character(Harvested))) %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  pull(Site) %>% 
  unique()

weather <- weather %>%
  filter(Site %in% ph_sites)
yield <- yield %>%
  filter(Site %in% ph_sites)
meta <- meta %>%
  filter(Site %in% ph_sites)


# Examine the quality of the temperature data -----------------------------
for (i in unique(weather$Year)) {
  filter(weather, Year == i) %>%
    mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
    group_by(Environment, Date) %>%
    summarise(TMIN = min(TEMP, na.rm = TRUE), 
              TMAX = max(TEMP, na.rm = TRUE)) %>%
    ungroup() %>%
    gather(Measure, TEMP, TMIN, TMAX) %>%
    ggplot(., aes(x = Date, y = TEMP, group = Measure, colour = Measure)) + 
      theme_bw() + geom_line() + facet_wrap(~ Environment) +
      scale_colour_manual(values = c("TMIN" = "blue", "TMAX" = "red"))
  ggsave(paste0("figures/munge/weather_unimputed_TEMP_", i, ".pdf"), height = 8, 
         width = 10, units = "in", dpi = 300)
}

# 2016 data look suspicious: late season temperatures never decline
weather <- weather %>%
  mutate(TEMP = if_else(Year == 2016, as.numeric(NA), TEMP))


# Collect days to impute temperature --------------------------------------
# Days that require minimum temperature imputation
min_temp_days <- weather %>%
  filter(Time >= hms("00:00:00") & Time <= hms("06:00:00")) %>%
  group_by(Year, Environment, Month, Day) %>%
  summarise(TEMP = sum(is.na(TEMP))) %>%
  ungroup() %>%
  filter(TEMP > 0) %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  pull(Key)
  
# Days that require maximum temperature imputation
max_temp_days <- weather %>%
  filter(Time >= hms("12:00:00") & Time <= hms("18:00:00")) %>%
  group_by(Year, Environment, Month, Day) %>%
  summarise(TEMP = sum(is.na(TEMP))) %>%
  ungroup() %>%
  filter(TEMP > 0) %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  pull(Key)

# Days that are completely missing from the weather files
missing_temp <- weather %>%
  split(., .$Site) %>%
  map(function(df) {
    days <- with(df, paste(Year, Month, Day, sep = "-")) %>%
      ymd() %>% unique()
    temp <- filter(yield, Year == df$Year[1], Environment == df$Environment[1])
    first <- min(temp$Planted, na.rm = TRUE)
    last <- max(temp$Harvested, na.rm = TRUE)
    missing <- seq(first, last, 1)
    missing <- missing[!(missing %in% days)]
    if (length(missing) == 0) return(NULL)
    paste(df$Environment[1], str_replace_all(missing, "-", "_"), sep = "_")
  })
missing_temp <- missing_temp[!sapply(missing_temp, is.null)]
missing_temp <- unlist(missing_temp, use.names = FALSE)

# Days from 2017 (for which we have absolutely no weather data)
days_2017 <- yield %>%
  filter(Year == 2017) %>%
  group_by(Year, Environment) %>%
  summarise(Planted = min(Planted, na.rm = TRUE), 
            Harvested = max(Harvested, na.rm = TRUE)) %>%
  ungroup() %>%
  by_row(function(r) {
    seq(r$Planted[1], r$Harvested[1], 1)
  }, .to = "Date") %>%
  unnest(Date) %>%
  mutate(Key = paste(Environment, Date, sep = "_") %>%
           str_replace_all(., "-", "_")) %>%
  pull(Key)


# Reformat and download PRISM data ----------------------------------------
min_temp <- tibble(Key = c(min_temp_days, missing_temp, days_2017)) %>%
  separate(Key, c("Environment", "Year", "Month", "Day"), sep = "_", remove = FALSE) %>%
  mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
  mutate_at(c("Year", "Month", "Day"), as.integer) %>%
  inner_join(., meta, by = c("Year", "Environment"))
max_temp <- tibble(Key = c(max_temp_days, missing_temp, days_2017)) %>%
  separate(Key, c("Environment", "Year", "Month", "Day"), sep = "_", remove = FALSE) %>%
  mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
  mutate_at(c("Year", "Month", "Day"), as.integer) %>%
  inner_join(., meta, by = c("Year", "Environment"))

options(prism.path = "~/prismtmp")
get_prism_dailys(type = "tmin", dates = unique(min_temp$Date), keepZip = FALSE)
get_prism_dailys(type = "tmax", dates = unique(max_temp$Date), keepZip = FALSE)


# Get temperature from the nearest PRISM grid -----------------------------
# Compute great circle distance between two points on earth's surface
source("src/haversine.R")

# Throw-away function for renaming columns
foo <- function(x) "Value"

# Create a cluster
cl <- makeCluster(detectCores())
clusterEvalQ(cl, { library(tidyverse); library(purrrlyr); library(prism); 
  source("src/haversine.R"); foo <- function(x) "Value" })

# Daily minimum temperature
split_idx <- clusterSplit(cl, 1:nrow(min_temp))
split_idx <- rep(1:detectCores(), times = sapply(split_idx, length))

min_temp <- min_temp %>%
  mutate(Date = str_replace_all(Date, "-", "")) %>%
  split(., split_idx)

min_temp <- parLapply(cl, min_temp, function(df) {
  by_row(df, function(r) {
    grep(ls_prism_data()[, 1], pattern = r$Date, value = TRUE) %>%
      grep(., pattern = "tmin", value = TRUE) %>%
      prism_stack() %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= r$Longitude - 2, Longitude <= r$Longitude + 2,
             Latitude >= r$Latitude - 2, Latitude <= r$Latitude + 2) %>%
      mutate(Distance = haversine(Latitude, r$Latitude, 
                                  Longitude, r$Longitude)) %>%
      rename_if(str_detect(names(.), "_"), funs(foo)) %>%
      arrange(Distance) %>%
      slice(1L) %>%
      pull(Value)
  }, .collate = "row", .to = "TMIN")
}) %>%
  bind_rows()
write_rds(min_temp, "data/weather/tmin_prism.rds")

# Daily maximum temperature
split_idx <- clusterSplit(cl, 1:nrow(max_temp))
split_idx <- rep(1:detectCores(), times = sapply(split_idx, length))

max_temp <- max_temp %>%
  mutate(Date = str_replace_all(Date, "-", "")) %>%
  split(., split_idx)

max_temp <- parLapply(cl, max_temp, function(df) {
  by_row(df, function(r) {
    grep(ls_prism_data()[, 1], pattern = r$Date, value = TRUE) %>%
      grep(., pattern = "tmax", value = TRUE) %>%
      prism_stack() %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= r$Longitude - 2, Longitude <= r$Longitude + 2,
             Latitude >= r$Latitude - 2, Latitude <= r$Latitude + 2) %>%
      mutate(Distance = haversine(Latitude, r$Latitude, 
                                  Longitude, r$Longitude)) %>%
      rename_if(str_detect(names(.), "_"), funs(foo)) %>%
      arrange(Distance) %>%
      slice(1L) %>%
      pull(Value)
  }, .collate = "row", .to = "TMAX")
}) %>%
  bind_rows()
write_rds(max_temp, "data/weather/tmax_prism.rds")


# Collect days to impute precipitation ------------------------------------
ppt_days <- weather %>%
  group_by(Environment, Year, Month, Day) %>%
  summarise(PPT = sum(PPT, na.rm = TRUE),
            N = n()) %>%
  ungroup() %>%
  filter(N < 48) %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  pull(Key)

# Use days in `missing_temp` and `days_2017` as well


# Reformat and download PRISM data ----------------------------------------
total_ppt <- tibble(Key = c(ppt_days, missing_temp, days_2017)) %>%
  separate(Key, c("Environment", "Year", "Month", "Day"), sep = "_", remove = FALSE) %>%
  mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
  mutate_at(c("Year", "Month", "Day"), as.integer) %>%
  inner_join(., meta, by = c("Year", "Environment"))

get_prism_dailys(type = "ppt", dates = unique(total_ppt$Date), keepZip = FALSE)


# Get precipitation from the nearest PRISM grid ---------------------------
split_idx <- clusterSplit(cl, 1:nrow(total_ppt))
split_idx <- rep(1:detectCores(), times = sapply(split_idx, length))

total_ppt <- total_ppt %>%
  mutate(Date = str_replace_all(Date, "-", "")) %>%
  split(., split_idx)

total_ppt <- parLapply(cl, total_ppt, function(df) {
  by_row(df, function(r) {
    grep(ls_prism_data()[, 1], pattern = r$Date, value = TRUE) %>%
      grep(., pattern = "ppt", value = TRUE) %>%
      prism_stack() %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= r$Longitude - 2, Longitude <= r$Longitude + 2, 
             Latitude >= r$Latitude - 2, Latitude >= r$Latitude + 2) %>%
      mutate(Distance = haversine(Latitude, r$Latitude, 
                                  Longitude, r$Longitude)) %>%
      rename_if(str_detect(names(.), "_"), funs(foo)) %>%
      arrange(Distance) %>%
      slice(1L) %>%
      pull(Value)
  }, .collate = "row", .to = "PPT")
}) %>%
  bind_rows()
write_rds(total_ppt, "data/weather/ppt_prism.rds")

# Stop the cluster
stopCluster(cl)
rm(cl); gc()


# Combine (non-)imputed data ----------------------------------------------
# Summarize daily values
daily_tmin <- weather %>%
  group_by(Environment, Year, Month, Day) %>%
  summarise(TMIN = min(TEMP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  filter(!(Key %in% min_temp$Key)) %>%
  bind_rows(.,
            select(min_temp, -Date, -Latitude, -Longitude, -Site, -Kernels))

daily_tmax <- weather %>%
  group_by(Environment, Year, Month, Day) %>%
  summarise(TMAX = max(TEMP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  filter(!(Key %in% max_temp$Key)) %>%
  bind_rows(.,
            select(max_temp, -Date, -Latitude, -Longitude, -Site, -Kernels))

daily_ppt <- weather %>%
  group_by(Environment, Year, Month, Day) %>%
  summarise(PPT = sum(PPT, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  filter(!(Key %in% total_ppt$Key)) %>%
  bind_rows(., 
            select(total_ppt, -Date, -Latitude, -Longitude, -Site, -Kernels))

# Some days are still missing for TMIN and TMAX
for (i in which(is.infinite(daily_tmin$TMIN))) {
  r <- daily_tmin[i, ]
  r$Latitude <- meta$Latitude[r$Environment == meta$Environment &
                                r$Year == meta$Year]
  r$Longitude <- meta$Longitude[r$Environment == meta$Environment &
                                  r$Year == meta$Year]
  get_prism_dailys(type = "tmin", 
                   dates = ymd(paste(r$Year, r$Month, r$Day, sep = "-")), 
                   keepZip = FALSE)
  
  daily_tmin$TMIN[i] <- grep(ls_prism_data()[, 1], 
                             pattern = ymd(paste(r$Year, r$Month, r$Day, sep = "-")) %>%
                               str_replace_all(., "-", ""), value = TRUE) %>%
    grep(., pattern = "tmin", value = TRUE) %>%
    prism_stack() %>%
    raster::rasterToPoints() %>%
    as_tibble() %>%
    rename(Longitude = x, Latitude = y) %>%
    filter(Longitude >= r$Longitude - 2, Longitude <= r$Longitude + 2, 
           Latitude >= r$Latitude - 2, Latitude >= r$Latitude + 2) %>%
    mutate(Distance = haversine(Latitude, r$Latitude, 
                                Longitude, r$Longitude)) %>%
    rename_if(str_detect(names(.), "_"), funs(foo)) %>%
    arrange(Distance) %>%
    slice(1L) %>%
    pull(Value)
}

for (i in which(is.infinite(daily_tmax$TMAX))) {
  r <- daily_tmax[i, ]
  r$Latitude <- meta$Latitude[r$Environment == meta$Environment &
                                r$Year == meta$Year]
  r$Longitude <- meta$Longitude[r$Environment == meta$Environment &
                                  r$Year == meta$Year]
  get_prism_dailys(type = "tmax", 
                   dates = ymd(paste(r$Year, r$Month, r$Day, sep = "-")), 
                   keepZip = FALSE)
  
  daily_tmax$TMAX[i] <- grep(ls_prism_data()[, 1], 
                             pattern = ymd(paste(r$Year, r$Month, r$Day, sep = "-")) %>%
                               str_replace_all(., "-", ""), value = TRUE) %>%
    grep(., pattern = "tmax", value = TRUE) %>%
    prism_stack() %>%
    raster::rasterToPoints() %>%
    as_tibble() %>%
    rename(Longitude = x, Latitude = y) %>%
    filter(Longitude >= r$Longitude - 2, Longitude <= r$Longitude + 2, 
           Latitude >= r$Latitude - 2, Latitude >= r$Latitude + 2) %>%
    mutate(Distance = haversine(Latitude, r$Latitude, 
                                Longitude, r$Longitude)) %>%
    rename_if(str_detect(names(.), "_"), funs(foo)) %>%
    arrange(Distance) %>%
    slice(1L) %>%
    pull(Value)
}

# Combine all variables
daily_tmin <- select(daily_tmin, -Key)
daily_tmax <- select(daily_tmax, -Key)
daily_ppt <- select(daily_ppt, -Key)

weather_imputed <- daily_tmin %>%
  inner_join(., daily_tmax, by = c("Environment", "Year", "Month", "Day")) %>%
  inner_join(., daily_ppt, by = c("Environment", "Year", "Month", "Day")) %>%
  select(Year, Environment, Month, Day, TMIN, TMAX, PPT) %>%
  arrange(Year, Environment, Month, Day)

write_rds(weather_imputed, "data/weather/weather_prism.rds")


# Re-examine the imputed temperatures -------------------------------------
for (i in unique(weather_imputed$Year)) {
  filter(weather_imputed, Year == i) %>%
    mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
    gather(Measure, TEMP, TMIN, TMAX) %>%
    ggplot(., aes(x = Date, y = TEMP, group = Measure, colour = Measure)) + 
      theme_bw() + geom_line() + facet_wrap(~ Environment) +
      scale_colour_manual(values = c("TMIN" = "blue", "TMAX" = "red"))
  ggsave(paste0("figures/munge/weather_imputed_TEMP_", i, ".pdf"), height = 8, 
         width = 10, units = "in", dpi = 300)
}
