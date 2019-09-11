library(tidyverse)
library(lubridate)
library(prism)
library(parallel)


# Download PRISM data -----------------------------------------------------
options(prism.path = "~/prismtmp")
get_prism_dailys("tmin", minDate = "1981-01-01", maxDate = "2018-12-31", keepZip = FALSE)
get_prism_dailys("tmax", minDate = "1981-01-01", maxDate = "2018-12-31", keepZip = FALSE)
get_prism_dailys("ppt", minDate = "1981-01-01", maxDate = "2018-12-31", keepZip = FALSE)


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


# Get site elevation ------------------------------------------------------
sites <- with(meta, paste(Environment, Year, sep = "_")) %>% unique() %>% sort()
prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
meta <- meta %>%
  filter(paste(Environment, Year, sep = "_") %in% sites) %>%
  rename(x = Longitude, y = Latitude) %>%
  select(x, y, Year, Environment)
elev <- elevatr::get_elev_point(as.data.frame(meta), prj = prj, src = "epqs")
meta <- meta %>%
  select(Year, Environment, x, y) %>%
  rename(Longitude = x, Latitude = y) %>%
  inner_join(., elev@data, by = c("Year", "Environment")) %>%
  select(-elev_units) %>%
  rename(Elevation = elevation)


# Create weather data files for LARS --------------------------------------
source("src/haversine.R")
foo <- function(x) "Value"
files <- ls_prism_data()[[1]]

cl <- makeCluster(detectCores())
clusterEvalQ(cl, { library(tidyverse); library(prism); options(prism.path = "~/prismtmp") })
clusterExport(cl, list("haversine", "foo", "files"))

sp <- split(meta, 1:nrow(meta))

ppt <- parLapply(cl, sp, function(df) {
  map_dbl(grep("ppt", files, value = TRUE), function(f) {
    prism_stack(f) %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= df$Longitude - 2, Longitude <= df$Longitude + 2,
             Latitude >= df$Latitude - 2, Latitude <= df$Latitude + 2) %>%
      mutate(Distance = haversine(Latitude, df$Latitude,
                                  Longitude, df$Longitude)) %>%
      rename_if(str_detect(names(.), "_"), funs(foo)) %>%
      arrange(Distance) %>%
      slice(1L) %>%
      pull(Value)
  })
})
names(ppt) <- sites
write_rds(ppt, "data/lars/ppt_raw.rds")

tmin <- parLapply(cl, sp, function(df) {
  map_dbl(grep("tmin", files, value = TRUE), function(f) {
    prism_stack(f) %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= df$Longitude - 2, Longitude <= df$Longitude + 2, 
             Latitude >= df$Latitude - 2, Latitude <= df$Latitude + 2) %>%
      mutate(Distance = haversine(Latitude, df$Latitude, 
                                  Longitude, df$Longitude)) %>%
      rename_if(str_detect(names(.), "_"), funs(foo)) %>%
      arrange(Distance) %>%
      slice(1L) %>%
      pull(Value)
  })
})
names(tmin) <- sites
write_rds(tmin, "data/lars/tmin_raw.rds")

tmax <- parLapply(cl, sp, function(df) {
  map_dbl(grep("tmax", files, value = TRUE), function(f) {
    prism_stack(f) %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= df$Longitude - 2, Longitude <= df$Longitude + 2, 
             Latitude >= df$Latitude - 2, Latitude <= df$Latitude + 2) %>%
      mutate(Distance = haversine(Latitude, df$Latitude, 
                                  Longitude, df$Longitude)) %>%
      rename_if(str_detect(names(.), "_"), funs(foo)) %>%
      arrange(Distance) %>%
      slice(1L) %>%
      pull(Value)
  })
})
names(tmax) <- sites
write_rds(tmax, "data/lars/tmax_raw.rds")

stopCluster(cl)


# Reformat the raw data ---------------------------------------------------
dates <- grep("ppt", files, value = TRUE) %>%
  str_split(., "_") %>% sapply(., function(x) x[5]) %>%
  ymd(.)
years <- year(dates)
jdays <- yday(dates)

for (i in seq_along(sites)) {
  tibble(Year = years, Day = jdays, TMIN = tmin[[i]], TMAX = tmax[[i]], 
         PPT = ppt[[i]]) %>%
    write_tsv(., paste0("data/lars/historical/", sites[i], "_historical.sr"), 
              col_names = FALSE)
}


# Construct site files ----------------------------------------------------
for (i in seq_along(sites)) {
  fname <- paste0("data/lars/historical/", sites[i], "_historical.st")
  cat("[SITE]", 
      sites[i], 
      "[LAT, LON and ALT]", 
      paste(meta$Latitude[i], meta$Longitude[i], meta$Elevation[i], collapse = "\t"), 
      "[WEATHER FILES]", 
      paste0("D:\\Aaron\\lars\\historical\\", sites[i], "_historical.sr"), 
      "[FORMAT]", 
      "YEAR JDAY MIN MAX RAIN RAD", 
      "[END]", file = fname, sep = "\n")
}
