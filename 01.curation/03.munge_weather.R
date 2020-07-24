### Download and process weather data from G2F repositories

library(tidyverse)
library(lubridate)


# Download files ----------------------------------------------------------
### 2017 data is currently downloaded from a private cooperator link

### See 01.munge_metadata.R for information about configuring icommands to
### download data from CyVerse

### This line will need to change based on where you installed icommands.
Sys.setenv("PATH" = paste("/mnt/01/amkusmec/bin/icommands", Sys.getenv("PATH"), sep = ":"))

setwd("data/weather")
system("iinit")
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/c._2014_weather_data/", 
              "_g2f_2014_weather_data_description.txt"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/c._2014_weather_data/", 
              "g2f_2014_weather_clean.csv"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Mar_2017/c._2015_weather_data/", 
              "_g2f_2015_weather_data_description.txt"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Mar_2017/c._2015_weather_data/", 
              "g2f_2015_weather_clean.csv"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "GenomesToFields_G2F_2016_Data_Mar_2018/c._2016_weather_data/", 
              "_g2f_2016_weather_data_description.txt"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "GenomesToFields_G2F_2016_Data_Mar_2018/c._2016_weather_data/", 
              "g2f_2016_weather_clean.csv"))
system("iexit")
setwd("../..")


# Load metadata -----------------------------------------------------------
# We need to look up time zones for 2016 sites for conversion of local datetime
# to UTC. We use the "fast" method from lutz because the accurate method requires
# GDAL>=2.0 which isn't available for Ubuntu 14.04.
meta <- read_csv("data/metadata/metadata_clean.csv") %>%
  mutate(Zone = lutz::tz_lookup_coords(Latitude, Longitude, warn = FALSE))

# Load and process files --------------------------------------------------
# Only keep weather data from sites in the metadata file
weather_2014 <- read_csv("data/weather/g2f_2014_weather_clean.csv") %>%
  select(`Experiment(s)`:`Dew point [C]`, `Solar Radiation [W/m2]`, 
         `Rainfall [mm]`) %>%
  rename(Environment = `Experiment(s)`, Day = `Day [Local]`,
         Month = `Month [Local]`, Year = `Year [Local]`, 
         Time = `Time [Local]`, UTC = `Datetime [UTC]`, 
         TEMP = `Temperature [C]`, DEW = `Dew point [C]`, 
         SR = `Solar Radiation [W/m2]`, PPT = `Rainfall [mm]`) %>%
  mutate(UTC = mdy_hm(UTC), 
         Environment = str_split(Environment, " ")) %>%
  select(-`Day of Year [Local]`) %>%
  unnest(Environment) %>%
  filter(str_detect(Environment, "[A-Z]{2}H")) %>%
  select(Environment, Day:PPT) %>%
  filter(Environment %in% meta$Environment, Year %in% meta$Year)

weather_2015 <- read_csv("data/weather/g2f_2015_weather_clean.csv") %>%
  select(`Experiment(s)`:`Dew Point [C]`, `Solar Radiation [W/m2]`, 
         `Rainfall [mm]`) %>%
  rename(Environment = `Experiment(s)`, Time = `Time [Local]`, 
         UTC = `Datetime [UTC]`, TEMP = `Temperature [C]`, 
         DEW = `Dew Point [C]`, SR = `Solar Radiation [W/m2]`, 
         PPT = `Rainfall [mm]`) %>%
  mutate(UTC = mdy_hm(UTC), 
         Environment = str_split(Environment, " ")) %>%
  select(-`Day of Year`) %>%
  unnest(Environment) %>%
  filter(str_detect(Environment, "[A-Z]{2}H")) %>%
  select(Environment, Day:PPT) %>%
  filter(Environment %in% meta$Environment, Year %in% meta$Year)

# The README for 2016 says that UTC datetime is in the file but this is not the
# case. We will need to convert local datetime to UTC using `force_tzs` for later
# standardization of SR values with observations from POWER.
weather_2016 <- read_csv("data/weather/g2f_2016_weather_clean.csv") %>%
  select(`Experiment(s)`:Year, `Time [Local]`:`Dew Point [C]`, 
         `Solar Radiation [W/m2]`, `Rainfall [mm]`) %>%
  rename(Environment = `Experiment(s)`, Time = `Time [Local]`, 
         TEMP = `Temperature [C]`, DEW = `Dew Point [C]`, 
         SR = `Solar Radiation [W/m2]`, PPT = `Rainfall [mm]`) %>%
  mutate(Environment = str_split(Environment, " "), 
         UTC = paste(Month, Day, Year, sep = "/") %>%
           paste(., Time) %>% mdy_hms(.)) %>%
  unnest(Environment) %>%
  filter(str_detect(Environment, "[A-Z]{2}H")) %>%
  inner_join(., meta, by = c("Year", "Environment")) %>%
  mutate(UTC = force_tzs(UTC, tzones = Zone, tzone_out = "UTC")) %>%
  select(Environment, Day:Time, UTC, TEMP:PPT)

### N.B. Weather data for 2017 is currently unavailable (1/18/2019).

### Units for SR are W/m^2 from G2F weather stations. Conversion to MJ/d/m^2 for
### evapotranspiration estimation comes from Figure 4.3.1.1 of
### https://power.larc.nasa.gov/new/files/POWER_Data_v8_methodology.pdf.
### kWh/d/m^2 = 3.6   MJ/d/m^2
### W/m^2     = 41.67 kWh/d/m^2
weather <- bind_rows(weather_2014, weather_2015, weather_2016) %>%
  mutate(SR = 3.6*(1/41.67)*SR)

write_rds(weather, "data/weather/weather_munged.rds")
