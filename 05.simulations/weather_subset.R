library(tidyverse)
library(lubridate)
library(prism)

options(prism.path = "~/prismtmp")
files <- ls_prism_data()[[1]]

long_min <- -94; long_max <- -93
lat_min <- 42; lat_max <- 43

tmin_files <- files[str_detect(files, "tmin")]
tmin <- map_df(tmin_files, function(f) {
    cat(f, "\r")
    prism_stack(f) %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= long_min, Longitude <= long_max, 
             Latitude >= lat_min, Latitude <= lat_max) %>%
      mutate(Date = str_extract(names(.), "[0-9]{8}")[3]) %>%
      rename_if(str_detect(names(.), "_"), function(x) "TMIN")
  })

tmax_files <- files[str_detect(files, "tmax")]
tmax <- map_df(tmax_files, function(f) {
    cat(f, "\r")
    prism_stack(f) %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= long_min, Longitude <= long_max, 
             Latitude >= lat_min, Latitude <= lat_max) %>%
      mutate(Date = str_extract(names(.), "[0-9]{8}")[3]) %>%
      rename_if(str_detect(names(.), "_"), function(x) "TMAX")
  })

ppt_files <- files[str_detect(files, "ppt")]
ppt <- map_df(ppt_files, function(f) {
    cat(f, "\r")
    prism_stack(f) %>%
      raster::rasterToPoints() %>%
      as_tibble() %>%
      rename(Longitude = x, Latitude = y) %>%
      filter(Longitude >= long_min, Longitude <= long_max, 
             Latitude >= lat_min, Latitude <= lat_max) %>%
      mutate(Date = str_extract(names(.), "[0-9]{8}")[3]) %>%
      rename_if(str_detect(names(.), "_"), function(x) "PPT")
  })

write_rds(list(tmin, tmax, ppt), "05.simulations/temp.rds")
