library(tidyverse)
library(readxl)


# Download the files ------------------------------------------------------
# 2017 data is currently downloaded through a private, cooperator link

### To download data follow these steps:
### 1) Download and install icommands. Information can be found at
###    https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands
###    For a Ubuntu server, run the following at the command prompt:
###      apt-get install https://files.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-icommands-4.1.10-ubuntu14-x86_64.deb
### 2) After installation, run `iinit` at the command prompt. Provide the
###    following responses to the prompts:
###      Enter the host name (DNS) of the server to connect to:data.cyverse.org
###      Enter the port number:1247
###      Enter your irods user name:anonymous

### This line will need to change based on where you installed icommands.
Sys.setenv("PATH" = paste("/mnt/01/amkusmec/bin/icommands", Sys.getenv("PATH"), sep = ":"))

setwd("data/metadata")
system("iinit")
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/z._2014_supplemental_info/", 
              "g2f_2014_field_characteristics.csv"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Mar_2017/z._2015_supplemental_info/", 
              "g2f_2015_field_metadata.csv"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "GenomesToFields_G2F_2016_Data_Mar_2018/z._2016_supplemental_info/", 
              "g2f_2016_field_metadata.csv"))
system("iexit")
setwd("../..")


# Process files -----------------------------------------------------------
meta_2014 <- read_csv("data/metadata/g2f_2014_field_characteristics.csv") %>%
  filter(Type == "hybrid") %>%
  rename(Environment = Experiment, Latitude = lat, Longitude = long, 
         Kernels = `Number kernels planted`) %>%
  mutate(Year = 2014,
         Environment = if_else(Environment == "G2FWI-HYB", "WIH1", Environment)) %>%
  select(Year, Environment, Latitude, Longitude, Kernels)

meta_2015 <- read_csv("data/metadata/g2f_2015_field_metadata.csv") %>%
  filter(Type == "Hybrid") %>%
  mutate(Year = 2015) %>%
  rename(Environment = Experiment, Latitude = `WS Lat`, Longitude = `WS Lon`, 
         Kernels = KernelsPerPlot) %>%
  select(Year, Environment, Latitude, Longitude, Kernels)

meta_2016 <- read_csv("data/metadata/g2f_2016_field_metadata.csv") %>%
  rename(Environment = `Experiment Code`,
         Latitude = `Weather station latitude (in decimal numbers NOT DMS)`,
         Longitude = `Weather station longitude (in decimal numbers NOT DMS)`, 
         Kernels = `Number kernels planted per plot (>200 seed/pack for cone planters)`) %>%
  mutate(Environment = str_replace(Environment, "--5/25 not done", ""),
         Year = 2016) %>%
  select(Year, Environment, Latitude, Longitude, Kernels)

# Includes two sites from NE that were sent via email
meta_2017 <- read_xlsx("data/metadata/metadata_2017_2018.xlsx") %>%
  rename(Year = Season, Environment = `Field-Location`, 
         Latitude = `Weather Station \r\nLatitude`, 
         Longitude = `Weather Station \r\nLogitude`) %>%
  select(Year, Environment, Latitude, Longitude) %>%
  filter(Year == 2017) %>%
  mutate(Kernels = NA) %>%
  bind_rows(., 
            tibble(Year = 2017, 
                   Environment = c("NEH3", "NEH4"), 
                   Latitude = c(40.9315, 40.9412), 
                   Longitude = c(-101.7670, -101.7659), 
                   Kernels = NA))

metadata <- bind_rows(meta_2014, meta_2015, meta_2016, meta_2017) %>%
  filter(!is.na(Latitude), !is.na(Longitude), !str_detect(Environment, "ON"))

write_csv(metadata, "data/metadata/metadata_clean.csv")


sites <- read_rds("data/phenotype/yield_agron0.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015") %>%
  pull(Site) %>% unique()
trials <- metadata %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  filter(Site %in% sites) %>%
  mutate(Latitude = 0.5*(Latitude %/% 0.5), 
         Longitude = 0.5*(Longitude %/% 0.5)) %>%
  count(Latitude, Longitude) %>%
  mutate(n = factor(n))
states <- map_data("state")
ggplot(states, aes(x = long, y = lat)) + theme_bw() + 
  geom_polygon(aes(group = group), fill = "white", colour = "black") + 
  geom_point(aes(x = Longitude, y = Latitude, size = n, fill = n), data = trials, 
             shape = 21, colour = "black") + 
  scale_fill_brewer(type = "seq", palette = "YlOrRd") +
  labs(x = "", y = "", fill = "# Trials", size = "# Trials")
ggsave("figures/trial_map.pdf", width = 7, height = 4, units = "in", dpi = 300)
