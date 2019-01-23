library(tidyverse)
library(purrrlyr)
library(lubridate)
library(jsonlite)


# Keep sites for which we have complete information -----------------------
weather <- read_rds("data/weather/weather_munged.rds") %>%
  select(-TEMP, -PPT, -DEW)
yield <- read_rds("data/phenotype/yield_munged.rds")
meta <- read_csv("data/metadata/metadata_clean.csv")

weather_sites <- with(weather, paste(Environment, Year, sep = "_")) %>% unique()
yield_sites <- with(yield, paste(Environment, Year, sep = "_")) %>% unique()
meta_sites <- with(meta, paste(Environment, Year, sep = "_")) %>% unique()
sites <- intersect(weather_sites, intersect(yield_sites, meta_sites))

weather <- weather %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites)
yield <- yield %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites)
meta <- meta %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites)


# Examine the quality of the solar radiation data -------------------------
# Based on conversion to POWER-reported values where the POWER values are the 
# average of seven daily values taken at 3, 6, 9, 12, 15, 18, and 21 hours GMT.
for (i in unique(weather$Year)) {
  filter(weather, Year == i) %>%
    filter(hour(UTC) %in% seq(3, 21, 3)) %>%
    group_by(Environment, Month, Day) %>%
    summarise(SR = mean(SR, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Date = ymd(paste(i, Month, Day, sep = "-"))) %>%
    ggplot(., aes(x = Date, y = SR)) + theme_bw() +
      geom_line() + facet_wrap(~ Environment)
  ggsave(paste0("figures/munge/weather_unimputed_SR_", i, ".pdf"), height = 8, 
         width = 10, units = "in", dpi = 300)
}

# SR data appear to have many zeros, completely missing days, and strange spikes.
# It's probably better to just gather data from POWER exclusively.


# Prepare the dates and locations -----------------------------------------
sr <- read_rds("data/weather/weather_prism.rds") %>%
  select(Year:Day) %>%
  mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
  group_by(Year, Environment) %>%
  summarise(Planted = min(Date), 
            Harvested = max(Date)) %>%
  ungroup() %>%
  inner_join(., meta, by = c("Year", "Environment"))


# Retrieve data from POWER ------------------------------------------------
# Operates through the POWER API (version 1.1.0; 23 January 2019)
# Retrieves data measured by satellite and interpolated on a daily 0.5 x 0.5
# degree grid.
# 
# Variables gathered include:
#  - ALLSKY_SFC_SW_DWN = all sky insolation incident on a horizontal surface (MJ/m^2/d)
#  - CLRSKY_SFC_SW_DWN = clear sky insolation incident on a horizontal surface (MJ/m^2/d)
#  - T2M_MAX = maximum temperature at 2 m (C)
#  - T2M_MIN = minimum temperature at 2 m (C)
#  - T2MDEW = dew/frost point at 2 m (C)
#  - WS2M = wind speed at 2 m (m/s)

# Function to construct the URL for data retrieval
POWER <- function(start_date, end_date, lat, lon) {
  paste0("https://power.larc.nasa.gov/cgi-bin/v1/DataAccess.py?",
         "request=execute&identifier=SinglePoint",
         "&parameters=ALLSKY_SFC_SW_DWN,CLRSKY_SFC_SW_DWN,T2M_MAX,T2M_MIN,T2MDEW,WS2M",
         "&startDate=", start_date,
         "&endDate=", end_date,
         "&lat=", lat,
         "&lon=", lon,
         "&userCommunity=AG",
         "&tempAverage=DAILY",
         "&outputList=JSON")
}


# Retrieve POWER data for all sites ---------------------------------------
sr <- sr %>%
  by_row(function(r) {
    temp <- POWER(start_date = str_replace_all(r$Planted[1], "-", ""), 
                  end_date = str_replace_all(r$Harvested[1], "-", ""), 
                  lat = r$Latitude, lon = r$Longitude) %>% 
      fromJSON()
    do.call("rbind", temp$features$properties$parameter) %>%
      t() %>%
      as_tibble(rownames = "Date") %>%
      mutate(Date = ymd(Date), 
             Month = month(Date), 
             Day = day(Date)) %>%
      mutate_at(vars(-Date, -Month, -Day), 
                function(x) if_else(x == -99, as.numeric(NA), x))
  }) %>%
  unnest(.out) %>%
  select(Year, Environment, Month, Day, ALLSKY_SFC_SW_DWN:WS2M)

write_rds(sr, "data/weather/weather_power.rds")


# Examine the radiation values from POWER ---------------------------------
for (i in unique(sr$Year)) {
  filter(sr, Year == i) %>%
    mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
    ggplot(., aes(x = Date, y = ALLSKY_SFC_SW_DWN)) + 
      theme_bw() + geom_line() + facet_wrap(~ Environment)
  ggsave(paste0("figures/munge/weather_imputed_SR_", i, ".pdf"), height = 8, 
         width = 10, units = "in", dpi = 300)
}
