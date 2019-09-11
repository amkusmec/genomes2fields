library(tidyverse)
library(purrrlyr)
library(lubridate)
library(jsonlite)


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


# Functions to retrieve data from POWER -----------------------------------
# Operates through the POWER API (version 1.1.0; 23 January 2019)
# Retrieves data measured by satellite and interpolated on a daily 0.5 x 0.5
# degree grid.
#
# Variables gathered include:
#  - ALLSKY_SFC_SW_DWN = all sky insolation incident on a horizontal surface (MJ/m^2/d)

POWER <- function(start_date, end_date, lat, lon) {
  paste0("https://power.larc.nasa.gov/cgi-bin/v1/DataAccess.py?", 
         "request=execute&identifier=SinglePoint",
         "&parameters=ALLSKY_SFC_SW_DWN", 
         "&startDate=", start_date, 
         "&endDate=", end_date, 
         "&lat=", lat, 
         "&lon=", lon,
         "&userCommunity=AG", 
         "&tempAverage=DAILY", 
         "&outputList=JSON")
}


# Retrieve data -----------------------------------------------------------
sr <- meta %>%
  by_row(function(r) {
    temp <- POWER(start_date = "19840101", end_date = "20181231", 
                  lat = r$Latitude[1], lon = r$Longitude[1]) %>%
      fromJSON()
    t(temp$features$properties$parameter$ALLSKY_SFC_SW_DWN) %>%
      as_tibble(rownames = "Date") %>%
      mutate(Date = ymd(Date)) %>%
      mutate_at(vars(-Date), 
                function(x) if_else(x == -99, as.numeric(NA), x)) %>%
      rename(ALLSKY_SFC_SW_DWN = `1`)
  }) %>%
  unnest(.out) %>%
  mutate(Year = year(Date), 
         JDAY = yday(Date)) %>%
  select(Site, Year, JDAY, ALLSKY_SFC_SW_DWN) %>%
  arrange(Year, JDAY) %>%
  split(., .$Site)


# Combine with data from PRISM --------------------------------------------
sites <- names(sr)
for (s in sites) {
  read_tsv(paste0("data/lars/historical/", s, "_historical.sr"), col_names = FALSE) %>%
    inner_join(., sr[[s]], by = c("X1" = "Year", "X2" = "JDAY")) %>%
    select(-Site) %>%
    write_tsv(., paste0("data/lars/historical/", s, "_historical.sr"), 
              col_names = FALSE)
}

# Replace missing values
for (f in list.files("data/lars/historical/", "*\\.sr", full.names = TRUE)) {
  read_tsv(f, col_names = FALSE) %>%
    mutate(X3 = if_else(is.na(X3), -99.0, X3), 
           X4 = if_else(is.na(X4), -99.0, X4), 
           X5 = if_else(is.na(X5), -99.0, X5), 
           X6 = if_else(is.na(X6), -99.0, X6)) %>%
    write_tsv(., f, col_names = FALSE)
}
